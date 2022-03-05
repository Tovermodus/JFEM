package multigrid;

import linalg.*;

public class JacobiSmoother
	implements Smoother
{
	private final int runs;
	SparseMatrix Dinv;
	final double omega;
	
	public JacobiSmoother(final int runs, final double omega, final Decomposable A)
	{
		this.omega = omega;
		this.runs = runs;
		final DenseVector AD =
			A.diag();
		Dinv = new SparseMatrix(A.getShape());
		for (int i = 0; i < AD.size(); i++)
		{
			if (Math.abs(AD.at(i)) < 1e-10)
				AD.set(1, i);
			Dinv.set(1. / AD.at(i), i, i);
		}
	}
	
	@Override
	public Vector smooth(final VectorMultiplyable operator, final Vector rhs, Vector iterate,
	                     final boolean verbose, final String prefix)
	{
		if (verbose)
		{
			final Vector residual = rhs.sub(operator.mvMul(iterate));
			System.out.println(prefix + residual.euclidianNorm());
		}
		for (int i = 0; i < runs; i++)
		{
			final Vector residual = rhs.sub(operator.mvMul(iterate));
//			if (verbose)
//				System.out.println(prefix + residual.euclidianNorm());
			iterate = iterate.add(Dinv.mvMul(residual)
			                          .mul(omega));
		}
		if (verbose)
		{
			final Vector residual = rhs.sub(operator.mvMul(iterate));
			System.out.println(prefix + residual.euclidianNorm());
		}
		return iterate;
	}
}
