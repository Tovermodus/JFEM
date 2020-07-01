package multigrid;

import basic.FESpace;
import linalg.DoubleTensor;

public class GaussSeidelSmoother extends Smoother
{
	private DoubleTensor r;
	private DoubleTensor l;
	public GaussSeidelSmoother(FESpace g, String[] args)
	{
		super(g, args);
		l = g.getSystemMatrix().getDiagonalMatrix().add(g.getSystemMatrix().getLowerTriangleMatrix());
		r = g.getSystemMatrix().getUpperTriangleMatrix();
	}

	@Override
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		DoubleTensor ret;
		ret = getL().solve(rightHandSide.sub(getR().mvmul(iterate)));
		return ret;
	}
	
	public DoubleTensor getR()
	{
		return r;
	}
	
	public DoubleTensor getL()
	{
		return l;
	}
}
