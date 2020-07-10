package basic;

import linalg.CoordinateVector;

public class LagrangeNodeFunctional implements NodeFunctional
{
	private CoordinateVector point;
	public LagrangeNodeFunctional(CoordinateVector point)
	{
		this.point = point;
	}
	
	public CoordinateVector getPoint()
	{
		return point;
	}
	
	@Override
	public double evaluate(ScalarFunction func)
	{
		 return func.value(point);
	}
}
