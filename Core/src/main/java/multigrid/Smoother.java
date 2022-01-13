package multigrid;

import linalg.Vector;
import linalg.VectorMultiplyable;

public interface Smoother
{
	Vector smooth(VectorMultiplyable Operator, Vector rhs, Vector iterate);
}
