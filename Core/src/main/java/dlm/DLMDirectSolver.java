package dlm;

import linalg.BlockSparseMatrix;
import linalg.DenseVector;
import linalg.Vector;

import java.util.List;

public class DLMDirectSolver
	extends DLMSolver
{
	@Override
	protected Vector solve(final BlockSparseMatrix systemMatrix,
	                       final DenseVector rhs,
	                       final FluidIterate fluidState,
	                       final List<ParticleIterate> particleStates,
	                       final FluidSystem fluidSystem,
	                       final List<ParticleSystem> particleSystems, final double dt, final double t)
	{
		return systemMatrix.toSparse()
		                   .solve(rhs);
	}
}
