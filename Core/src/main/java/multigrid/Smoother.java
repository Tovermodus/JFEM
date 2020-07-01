package multigrid;

import basic.FESpace;
import linalg.DoubleTensor;

public class Smoother
{
	FESpace fESpace;
	String[] args;
	public Smoother(FESpace g, String[] args)
	{
		this.fESpace = g;
		this.args = args;
	}
	public static Smoother fromSmoother(Smoother s)
	{
		return new Smoother(s.fESpace,s.args);
	}
	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
	{
		throw new UnsupportedOperationException();
	}
}

