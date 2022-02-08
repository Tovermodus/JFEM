package multigrid;

import linalg.*;

public class ForwardBackwardGaussSeidelSmoother
	implements Smoother
{
	private final int runs;
	LowerTriangularSparseMatrix L;
	Matrix D;
	UpperTriangularSparseMatrix U;
	
	public ForwardBackwardGaussSeidelSmoother(final int runs, final Decomposable A)
	{
		
		this.runs = runs;
		L = new LowerTriangularSparseMatrix(A);
		D = new SparseMatrix(A.getDiagonalMatrix());
		U = new UpperTriangularSparseMatrix(A);
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
			iterate = iterate.add(L.solve(D.mvMul(U.solve(residual))));
		}
		if (verbose)
		{
			final Vector residual = rhs.sub(operator.mvMul(iterate));
			System.out.println(prefix + residual.euclidianNorm());
		}
		return iterate;
	}
}
