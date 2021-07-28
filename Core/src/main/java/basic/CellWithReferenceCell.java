package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface CellWithReferenceCell<CT extends CellWithReferenceCell<CT, FT, ET>, FT extends FaceWithReferenceFace<CT, FT, ET>,
	ET extends EdgeWithReferenceEdge<CT, FT, ET>> extends Cell<CT,FT,ET>
{
	CT getReferenceCell();
	
	CoordinateMatrix transformationGradientFromReferenceCell(CoordinateVector pos);
	CoordinateMatrix transformationGradientToReferenceCell(CoordinateVector pos);
	CoordinateVector transformToReferenceCell(CoordinateVector pos);
	CoordinateVector transformFromReferenceCell(CoordinateVector pos);
}
