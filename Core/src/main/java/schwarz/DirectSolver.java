package schwarz;

import linalg.DenseMatrix;
import linalg.Matrix;
import linalg.Vector;

public class DirectSolver
	implements SystemSolver<Matrix>
{
	@Override
	public Vector solve(final Matrix system, final Vector rhs)
	{
		final DenseMatrix localOp = new DenseMatrix(system);
		return localOp.solve(rhs);
	}
}
