package dlm;

import linalg.BlockSparseMatrix;
import linalg.DenseVector;
import linalg.Vector;

import java.util.List;

public abstract class DLMSolver
{
	protected abstract Vector solve(BlockSparseMatrix systemMatrix, DenseVector rhs, FluidIterate fluidState,
	                                List<ParticleIterate> particleStates, FluidSystem fluidSystem,
	                                List<ParticleSystem> particleSystems, double dt, double t);
}
