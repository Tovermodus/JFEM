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
		Linv = new DenseMatrix(A.getLowerTriangleMatrix()
		                        .add(A.getDiagonalMatrix())).inverse();
		D = new SparseMatrix(A.getDiagonalMatrix());
		R = A.getUpperTriangleMatrix();
	}
	
	@Override
	public Vector smooth(final VectorMultiplyable operator, final Vector rhs, Vector iterate)
	{
		for (int i = 0; i < runs; i++)
		{
			iterate = Linv.mvMul(rhs.sub(R.mvMul(iterate)));
		}
		return iterate;
	}
}
