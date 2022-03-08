package schwarz;

import linalg.Vector;
import linalg.VectorMultiplyable;

public interface SystemSolver<OT extends VectorMultiplyable>
{
	Vector solve(OT system, Vector rhs, int patch);
}
