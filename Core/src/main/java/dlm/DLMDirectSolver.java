package dlm;

import linalg.BlockSparseMatrix;
import linalg.DenseVector;

import java.util.List;

public class DLMDirectSolver
	extends DLMSolver
{
	@Override
	protected DenseVector solve(final BlockSparseMatrix systemMatrix,
	                            final DenseVector rhs,
	                            final FluidIterate fluidState,
	                            final List<ParticleIterate> particleStates,
	                            final FluidSystem fluidSystem,
	                            final List<ParticleSystem> particleSystems)
	{
		return systemMatrix.toSparse()
		                   .solve(rhs);
	}
}
