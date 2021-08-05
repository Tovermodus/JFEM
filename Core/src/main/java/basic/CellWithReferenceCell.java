package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface CellWithReferenceCell<CT extends CellWithReferenceCell<CT, FT>, FT extends FaceWithReferenceFace<CT, FT>
	> extends Cell<CT,FT>
{
	CT getReferenceCell();
	
	CoordinateMatrix transformationGradientFromReferenceCell(CoordinateVector pos);
	CoordinateMatrix transformationGradientToReferenceCell(CoordinateVector pos);
	CoordinateVector transformToReferenceCell(CoordinateVector pos);
	CoordinateVector transformFromReferenceCell(CoordinateVector pos);
}
