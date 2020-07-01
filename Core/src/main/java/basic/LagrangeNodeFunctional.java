package basic;

import linalg.DoubleTensor;

public class LagrangeNodeFunctional extends ScalarNodeFunctional
{
	private DoubleTensor point;
	public LagrangeNodeFunctional(DoubleTensor point)
	{
		this.point = point;
	}
	
	public DoubleTensor getPoint()
	{
		return point;
	}
	
	@Override
	public double evaluate(ScalarFunction func)
	{
		return func.value(point);
	}
}
