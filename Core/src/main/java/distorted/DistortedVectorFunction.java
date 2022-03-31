package distorted;

import basic.VectorFunction;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

public interface DistortedVectorFunction
	extends ReferenceCellFunction<DistortedCell, DistortedFace, CoordinateVector, CoordinateMatrix,
	CoordinateTensor>,
	VectorFunction
{
}
