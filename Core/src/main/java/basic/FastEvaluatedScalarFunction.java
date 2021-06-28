package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface FastEvaluatedScalarFunction extends ScalarFunction
{
	double fastValue(CoordinateVector pos);
	double[] fastGradient(CoordinateVector pos);
	double[][] fastHessian(CoordinateVector pos);
	
	@Override
	default Double value(CoordinateVector pos)
	{
		return fastValue(pos);
	}
	
	@Override
	default CoordinateVector gradient(CoordinateVector pos)
	{
		return new CoordinateVector(fastGradient(pos));
	}
	
	@Override
	default CoordinateMatrix hessian(CoordinateVector pos)
	{
		return new CoordinateMatrix(fastHessian(pos));
	}
}
