package mixed;

import basic.*;
import linalg.CoordinateVector;

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
	
	public MixedFunctionOnCells()
	{
		super();
	}
	
	@Override
	public MixedValue valueInCell(final CoordinateVector pos, final CT cell)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public MixedGradient gradientInCell(final CoordinateVector pos, final CT cell)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
