package dlm;

import linalg.BlockSparseMatrix;
import linalg.DenseVector;

import java.util.List;

public abstract class DLMSolver
{
	protected abstract DenseVector solve(BlockSparseMatrix systemMatrix, DenseVector rhs, FluidIterate fluidState,
	                                     List<ParticleIterate> particleStates, FluidSystem fluidSystem,
	                                     List<ParticleSystem> particleSystems);
}
