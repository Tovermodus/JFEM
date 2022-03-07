package schwarz;

import linalg.IterativeSolver;
import linalg.Matrix;
import linalg.Vector;

public class CGSolver
	implements SystemSolver<Matrix>
{
	double tol;
	
	public CGSolver(final double tol)
	{
		this.tol = tol;
	}
	
	@Override
	public Vector solve(final Matrix system, final Vector rhs)
	{
		final IterativeSolver it = new IterativeSolver(true);
		return it.solveCG(system, rhs, tol);
	}
}
