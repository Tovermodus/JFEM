package mixed;

import basic.*;

public abstract class MixedFunctionOnCells<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	extends MixedFunction
	implements FunctionOnCells<CT,
	FT, MixedValue,
	MixedGradient,
	MixedHessian>
{
	public MixedFunctionOnCells(final ScalarFunctionOnCells<CT, FT> scalarFunction,
	                            final VectorFunctionOnCells<CT, FT> vectorFunction)
	{
		super(scalarFunction, vectorFunction);
	}
}
