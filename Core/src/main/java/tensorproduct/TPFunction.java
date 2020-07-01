package tensorproduct;

import basic.ScalarFunction;
import linalg.DoubleTensor;

public class TPFunction extends ScalarFunction
{
	Function1D xFunction;
	Function1D yFunction;
	public TPFunction(Function1D xFunction, Function1D yFunction)
	{
		this.xFunction = xFunction;
		this.yFunction = yFunction;
	}
	@Override
	public double value(DoubleTensor pos)
	{
		return xFunction.value(pos.x())*yFunction.value(pos.y());
	}

	@Override
	public DoubleTensor derivative(DoubleTensor pos)
	{
		return DoubleTensor.vectorFromValues(xFunction.derivative(pos.x())*yFunction.value(pos.y()),
			xFunction.value(pos.x())*yFunction.derivative(pos.y()));
	}
}
