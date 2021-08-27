package basic;

import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface FastEvaluatedScalarFunction extends ScalarFunction
{
	double fastValue(CoordinateVector pos);
	
	double[] fastGradient(CoordinateVector pos);
	
	double[][] fastHessian(CoordinateVector pos);
	
	@Override
	default Double value(final CoordinateVector pos)
	{
		return fastValue(pos);
	}
	
	@Override
	default CoordinateVector gradient(final CoordinateVector pos)
	{
		return new CoordinateVector(fastGradient(pos));
	}
	
	@Override
	default CoordinateMatrix hessian(final CoordinateVector pos)
	{
		return new CoordinateDenseMatrix(fastHessian(pos));
	}
}
