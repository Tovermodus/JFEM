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
}
