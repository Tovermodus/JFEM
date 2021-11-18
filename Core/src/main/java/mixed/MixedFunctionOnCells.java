package mixed;

import basic.*;

public interface MixedFunctionOnCells<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	extends MixedFunction, FunctionOnCells<CT, FT, MixedValue, MixedGradient, MixedHessian>
{
	@Override
	ScalarFunctionOnCells<CT, FT> getPressureFunction();
	
	@Override
	VectorFunctionOnCells<CT, FT> getVelocityFunction();
}
