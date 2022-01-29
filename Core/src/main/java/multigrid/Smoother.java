package multigrid;

import linalg.Vector;
import linalg.VectorMultiplyable;

public interface Smoother
{
	default Vector smooth(final VectorMultiplyable Operator, final Vector rhs, final Vector iterate)
	{
		return smooth(Operator, rhs, iterate, false);
	}
	
	Vector smooth(VectorMultiplyable Operator, Vector rhs, Vector iterate, boolean verbose);
}
