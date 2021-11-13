package basic;

import linalg.CoordinateVector;

public interface Function<valueT, gradientT, hessianT>
{
	int getDomainDimension(); // should return -1 if any dimension is accepted
	
	valueT defaultValue();
	
	gradientT defaultGradient();
	
	hessianT defaultHessian();
	
	valueT value(CoordinateVector pos);
	
	default gradientT gradient(final CoordinateVector pos)
	{
		throw new UnsupportedOperationException();
	}
	
	default hessianT hessian(final CoordinateVector pos)
	{
		throw new UnsupportedOperationException();
	}
	
	default FunctionSignature getFunctionSignature()
	{
		return new FunctionSignature(defaultValue().getClass(), defaultGradient().getClass(),
		                             defaultHessian().getClass());
	}
	
	Function<valueT, gradientT, hessianT> concatenateWith(VectorFunction f);
	
	default <CT extends Cell<CT, FT>, FT extends Face<CT, FT>> FunctionOnCells<CT, FT, valueT, gradientT,
		hessianT> concatenateWithOnCells(
		final VectorFunctionOnCells<CT, FT> f)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (f.getRangeDimension() != getDomainDimension())
				throw new IllegalArgumentException("Inner function has the wrong range");
		}
		
		final Function<valueT, gradientT, hessianT> function = concatenateWith((VectorFunction) f);
		return new FunctionOnCells<CT, FT, valueT, gradientT, hessianT>()
		{
			@Override
			public valueT valueInCell(final CoordinateVector pos, final CT cell)
			{
				return function.value(pos);
			}
			
			@Override
			public gradientT gradientInCell(final CoordinateVector pos, final CT cell)
			{
				return function.gradient(pos);
			}
			
			@Override
			public int getDomainDimension()
			{
				return function.getDomainDimension();
			}
			
			@Override
			public valueT defaultValue()
			{
				return function.defaultValue();
			}
			
			@Override
			public gradientT defaultGradient()
			{
				return function.defaultGradient();
			}
			
			@Override
			public hessianT defaultHessian()
			{
				return function.defaultHessian();
			}
			
			@Override
			public valueT value(final CoordinateVector pos)
			{
				return function.value(pos);
			}
			
			@Override
			public gradientT gradient(final CoordinateVector pos)
			{
				return function.gradient(pos);
			}
			
			@Override
			public Function<valueT, gradientT, hessianT> concatenateWith(final VectorFunction f)
			{
				return function.concatenateWith(f);
			}
		};
	}
}
