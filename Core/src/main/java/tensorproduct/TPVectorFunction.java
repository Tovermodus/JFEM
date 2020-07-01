package tensorproduct;

import basic.TensorFunction;
import linalg.DoubleTensor;

public class TPVectorFunction extends TensorFunction
{
	Function1D xFunctionx;
	Function1D xFunctiony;
	Function1D yFunctionx;
	Function1D yFunctiony;
	public TPVectorFunction(Function1D xFunctionx, Function1D xFunctiony, Function1D yFunctionx,
	                        Function1D yFunctiony)
	{
		this.xFunctionx = xFunctionx;
		this.xFunctiony = xFunctiony;
		this.yFunctionx = yFunctionx;
		this.yFunctiony = yFunctiony;
	}
	@Override
	public DoubleTensor value(DoubleTensor pos)
	{
		return DoubleTensor.vectorFromValues(xFunctionx.value(pos.x())*xFunctiony.value(pos.y()),
			yFunctionx.value(pos.x())*yFunctiony.value(pos.y()));
	}

	@Override
	public DoubleTensor derivative(DoubleTensor pos)
	{
		/*
		∂_xf_1 ∂_yf_1
		∂_xf_2 ∂_yf_2
		 */
		return DoubleTensor.squareMatrixFromValues(xFunctionx.derivative(pos.x())*xFunctiony.value(pos.y()),
			xFunctionx.value(pos.x())*xFunctiony.derivative(pos.y()),
			yFunctionx.derivative(pos.x())*yFunctiony.value(pos.y()),
			yFunctionx.value(pos.x())*yFunctiony.derivative(pos.y()));
	}
}
