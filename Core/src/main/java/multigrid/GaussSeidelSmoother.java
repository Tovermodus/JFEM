package multigrid;

import linalg.*;

public class GaussSeidelSmoother
	implements Smoother
{
	private final int runs;
	Matrix Linv;
	Matrix D;
	Matrix R;
	
	public GaussSeidelSmoother(final int runs, final Decomposable A)
	{
		
		this.runs = runs;
		Linv = new DenseMatrix(A.getStrictlyLowerTriangleMatrix()
		                        .add(A.getDiagonalMatrix())).inverse();
		D = new SparseMatrix(A.getDiagonalMatrix());
		R = A.getStrictlyUpperTriangleMatrix();
	}
	
	@Override
	public Vector smooth(final VectorMultiplyable operator, final Vector rhs, Vector iterate, final boolean verbose)
	{
		for (int i = 0; i < runs; i++)
		{
			iterate = Linv.mvMul(rhs.sub(R.mvMul(iterate)));
		}
		return iterate;
	}
}
