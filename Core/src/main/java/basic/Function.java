package basic;

import linalg.CoordinateVector;

public interface Function<valueT, gradientT, hessianT>
{
	int getDomainDimension(); // should return -1 if any dimension is accepted
	valueT defaultValue();
	gradientT defaultGradient();
	hessianT defaultHessian();
	valueT value(CoordinateVector pos);
	default gradientT gradient(CoordinateVector pos)
	{
		throw new UnsupportedOperationException();
	}
	default hessianT hessian(CoordinateVector pos)
	{
		throw new UnsupportedOperationException();
	}
	default FunctionIdentifier<?> createIdentifier(int number)
	{
		return new FunctionIdentifier<>(this.getClass(), number);
	}
}
