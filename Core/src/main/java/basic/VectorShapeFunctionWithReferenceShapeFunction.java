package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

public interface VectorShapeFunctionWithReferenceShapeFunction<CT extends CellWithReferenceCell<CT,FT>,
	FT extends FaceWithReferenceFace<CT,FT>>
	extends ShapeFunctionWithReferenceShapeFunction<CT,FT, CoordinateVector, CoordinateMatrix, CoordinateTensor>,
		VectorShapeFunction<CT,FT>
{
	@Override
	VectorShapeFunctionWithReferenceShapeFunction<CT, FT> createReferenceShapeFunctionRelativeTo(CT cell);
	
	@Override
	VectorShapeFunctionWithReferenceShapeFunction<CT, FT> createReferenceShapeFunctionRelativeTo(FT face);
	
}
