package distorted;

import basic.ScalarFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface DistortedScalarFunction
	extends DistortedFunction<Double, CoordinateVector, CoordinateMatrix>, ScalarFunction
{
}
