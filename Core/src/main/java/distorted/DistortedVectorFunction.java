package distorted;

import basic.VectorFunction;
import distorted.geometry.DistortedCell;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface DistortedVectorFunction extends VectorFunction
{
	
	CoordinateVector valueOnReferenceCell(final CoordinateVector pos, DistortedCell cell);
	
	CoordinateMatrix gradientOnReferenceCell(final CoordinateVector pos, DistortedCell cell);
}
