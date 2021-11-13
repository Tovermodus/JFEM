package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

public interface VectorFunctionOnCells<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	extends FunctionOnCells<CT,
	FT, CoordinateVector, CoordinateMatrix, CoordinateTensor>, VectorFunction
{
}
