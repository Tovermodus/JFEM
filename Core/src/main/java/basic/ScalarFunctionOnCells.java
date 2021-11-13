package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface ScalarFunctionOnCells<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	extends FunctionOnCells<CT,
	FT, Double, CoordinateVector, CoordinateMatrix>, ScalarFunction
{
}
