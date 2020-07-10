package basic;

import linalg.CoordinateVector;

public interface Function<valueT, gradientT, hessianT>
{
	int getDomainDimension(); // should return -1 if any dimension is accepted
	valueT value(CoordinateVector pos);
	default gradientT gradient(CoordinateVector pos)
	{
		throw new UnsupportedOperationException();
	}
	default hessianT hessian(CoordinateVector pos)
	{
		throw new UnsupportedOperationException();
	}
}
