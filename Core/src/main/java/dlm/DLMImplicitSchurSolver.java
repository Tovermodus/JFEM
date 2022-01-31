package dlm;

import linalg.BlockSparseMatrix;
import linalg.DenseVector;
import linalg.IterativeImplicitSchur;
import linalg.IterativeSolver;

import java.util.List;

public class DLMImplicitSchurSolver
	extends DLMSolver
{
	IterativeImplicitSchur schur;
	IterativeSolver it = new IterativeSolver(true);
	
	@Override
	protected DenseVector solve(final BlockSparseMatrix systemMatrix,
	                            final DenseVector rhs,
	                            final FluidIterate fluidState,
	                            final List<ParticleIterate> particleStates,
	                            final FluidSystem fluidSystem,
	                            final List<ParticleSystem> particleSystems)
	{
		
		it.showProgress = true;
		if (schur == null)
			schur = new IterativeImplicitSchur(systemMatrix);
		else
			schur.resetOffDiagonals(systemMatrix);
		return new DenseVector(it.solvePGMRES(systemMatrix, schur, rhs, 1e-7));
	}
}
