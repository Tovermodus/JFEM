package schwarz;

import linalg.Matrix;
import linalg.Vector;
import linalg.VectorMultiplyable;
import multigrid.Smoother;

public class SchwarzSmoother
	implements Smoother
{
	AbstractSchwarz<?, ?, Matrix> schwarz;
	int runs;
	
	public SchwarzSmoother(final int runs, final AbstractSchwarz<?, ?, Matrix> schwarz)
	{
		this.schwarz = schwarz;
		this.runs = runs;
	}
	
	@Override
	public Vector smooth(final VectorMultiplyable Operator,
	                     final Vector rhs,
	                     Vector iterate,
	                     final boolean verbose,
	                     final String prefix)
	{
		if (verbose)
			System.out.println(prefix + " " + schwarz.getGlobalOperator()
			                                         .mvMul(iterate)
			                                         .sub(rhs)
			                                         .euclidianNorm());
		for (int i = 0; i < runs; i++)
		{
			iterate = schwarz.getSubspaceCorrection()
			                 .apply(schwarz, iterate, rhs);
			if (verbose)
				System.out.println(prefix + " iter " + i + " " + schwarz.getGlobalOperator()
				                                                        .mvMul(iterate)
				                                                        .sub(rhs)
				                                                        .euclidianNorm());
		}
		return iterate;
	}
}
