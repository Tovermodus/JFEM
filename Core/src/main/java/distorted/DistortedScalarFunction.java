package distorted;

import basic.ScalarFunction;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface DistortedScalarFunction
	extends ReferenceCellFunction<DistortedCell, DistortedFace, Double, CoordinateVector, CoordinateMatrix>,
	ScalarFunction
{
}
