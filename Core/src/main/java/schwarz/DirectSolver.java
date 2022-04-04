package schwarz;

import linalg.DenseMatrix;
import linalg.Matrix;
import linalg.SparseMatrix;
import linalg.Vector;

import java.util.HashMap;
import java.util.Map;

public class DirectSolver
	implements SystemSolver<Matrix>
{
	Map<Integer, DenseMatrix> inverses;
	
	public DirectSolver()
	{
		inverses = new HashMap<>();
	}
	
	@Override
	public Vector solve(final Matrix system, final Vector rhs, final int patch)
	{
		if (!inverses.containsKey(patch))
		{
			final SparseMatrix s;
			if (system instanceof SparseMatrix)
				s = (SparseMatrix) system;
			else
				s = new SparseMatrix(system);
			final DenseMatrix inv = s.inverseNative();
			synchronized (this)
			{
				inverses.put(patch, inv);
			}
		}
		return inverses.get(patch)
		               .mvMul(rhs);
	}
}
