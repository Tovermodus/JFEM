package schwarz;

import linalg.DenseMatrix;
import linalg.Matrix;
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
			inverses.put(patch, new DenseMatrix(system).inverse());
		return inverses.get(patch)
		               .mvMul(rhs);
	}
}
