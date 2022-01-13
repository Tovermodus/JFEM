package multigrid;

import linalg.Vector;
import linalg.VectorMultiplyable;

public class RichardsonSmoother
	implements Smoother
{
	private final double omega;
	private final int runs;
	
	public RichardsonSmoother(final double omega, final int runs)
	{
		
		this.omega = omega;
		this.runs = runs;
	}
	
	@Override
	public Vector smooth(final VectorMultiplyable operator, final Vector rhs, Vector iterate)
	{
		for (int i = 0; i < runs; i++)
			iterate = iterate.add(rhs.sub(operator.mvMul(iterate))
			                         .mul(omega));
		return iterate;
	}
}
