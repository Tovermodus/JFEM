package multigrid;

import basic.MatrixFESpace;

public class GaussSeidelSmoother extends Smoother
{
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
