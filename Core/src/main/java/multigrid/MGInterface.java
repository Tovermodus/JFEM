package multigrid;

import linalg.Vector;
import linalg.VectorMultiplyable;

public interface MGInterface
{
	Vector mgStep(final int level, Vector guess, final Vector rhs);
	
	Vector vCycle(final Vector initialIterate, final Vector rhs);
	
	default Vector fullVCycleCorrection(final Vector rhs)
	{
		throw new UnsupportedOperationException("not yet implemented");
	}
	
	VectorMultiplyable getFinestSystem();
	
	Vector fullVCycleSolver(final Vector iterate);
}
