package multigrid;

import linalg.DenseMatrix;
import linalg.SparseMatrix;
import linalg.Vector;
import linalg.VectorMultiplyable;

public class GaussSeidelPreconditioner implements VectorMultiplyable
{
	DenseMatrix DU;
	DenseMatrix DL;
	SparseMatrix Dinv;
	
	public GaussSeidelPreconditioner(final SparseMatrix A)
	{
		System.out.println("crerating");
		DU =
			new DenseMatrix(A
				                .getDiagonalMatrix()
				                .add(SparseMatrix.identity(A.getCols()))
				                .add(A.getUpperTriangleMatrix())).inverse();
		System.out.println("crerating2");
		DL = new DenseMatrix(A
			                     .getDiagonalMatrix()
			                     .add(SparseMatrix.identity(A.getCols()))
			                     .add(A.getLowerTriangleMatrix())).inverse();
		System.out.println("crerating3");
		Dinv = new SparseMatrix(A.getRows(), A.getCols());
		for (int i = 0; i < A.getRows(); i++)
			Dinv.add((A.at(i, i) + 1), i, i);
	}
	
	@Override
	public int getVectorSize()
	{
		return DU.getCols();
	}
	
	@Override
	public int getTVectorSize()
	{
		return DU.getRows();
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		return DU.mvMul(Dinv.mvMul(DL.mvMul(vector)));
	}
	
	@Override
	public Vector tvMul(final Vector vector)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}

//	private final DoubleTensor r;
//	private final DoubleTensor l;
//	public GaussSeidelSmoother(MatrixFESpace<?,?,?> g, String[] args)
//	{
//		super(g, args);
//		l = g.getSystemMatrix().getDiagonalMatrix().add(g.getSystemMatrix().getLowerTriangleMatrix());
//		r = g.getSystemMatrix().getUpperTriangleMatrix();
//	}
//
//	@Override
//	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
//	{
//		DoubleTensor ret;
//		ret = getL().solve(rightHandSide.sub(getR().mvmul(iterate)));
//		return ret;
//	}
//
//	public DoubleTensor getR()
//	{
//		return r;
//	}
//
//	public DoubleTensor getL()
//	{
//		return l;
//	}
}
