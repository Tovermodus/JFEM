package dlm;

import linalg.*;

import java.util.List;

public class DLMImplicitSchurSolver
	extends DLMSolver
{
	IterativeImplicitSchur schur;
	IterativeSolver it = new IterativeSolver(true);
	
	@Override
	protected Vector solve(final BlockSparseMatrix systemMatrix,
	                       final DenseVector rhs,
	                       final FluidIterate fluidState,
	                       final List<ParticleIterate> particleStates,
	                       final FluidSystem fluidSystem,
	                       final List<ParticleSystem> particleSystems, final double dt, final double t)
	{
		
		it.showProgress = true;
		if (schur == null)
			schur = new IterativeImplicitSchur(systemMatrix);
		else
			schur.resetOffDiagonals(systemMatrix);
		return new DenseVector(it.solvePGMRES(systemMatrix, schur, rhs, 1e-7));
	}
}
