package distorted;

import basic.VectorFunction;
import basic.VectorFunctionOnCells;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

public interface DistortedVectorFunction
	extends DistortedFunction<CoordinateVector, CoordinateMatrix, CoordinateTensor>, VectorFunction
{
}
