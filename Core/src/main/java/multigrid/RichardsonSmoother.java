package multigrid;

import linalg.Vector;
import linalg.VectorMultiplyable;

public class RichardsonSmoother
	implements Smoother
{
	private final double omega;
	private final int runs;
	private final VectorMultiplyable preconditioner;
	
	public RichardsonSmoother(final double omega, final int runs)
	{
		
		this.omega = omega;
		this.runs = runs;
		preconditioner = null;
	}
	
	public RichardsonSmoother(final double omega, final int runs, final VectorMultiplyable preconditioner)
	{
		
		this.omega = omega;
		this.runs = runs;
		this.preconditioner = preconditioner;
	}
	
	@Override
	public Vector smooth(final VectorMultiplyable operator, final Vector rhs, Vector iterate,
	                     final boolean verbose, final String prefix)
	{
		for (int i = 0; i < runs; i++)
		{
			Vector residual = rhs.sub(operator.mvMul(iterate));
			if (verbose)
				System.out.println(prefix + residual.euclidianNorm());
			if (residual.euclidianNorm() < 1e-10)
			{
				if (verbose)
					System.out.println(prefix + "is smooth");
				break;
			}
			if (preconditioner != null)
				residual = preconditioner.mvMul(residual);
			iterate = iterate.add(residual.mul(omega));
		}
		return iterate;
	}
}
