package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface ScalarShapeFunctionWithReferenceShapeFunction<CT extends CellWithReferenceCell<CT,FT>,
	FT extends FaceWithReferenceFace<CT,FT>>
	extends ShapeFunctionWithReferenceShapeFunction<CT,FT, Double, CoordinateVector, CoordinateMatrix>,
		ScalarShapeFunction<CT,FT>
{
	@Override
	ScalarShapeFunctionWithReferenceShapeFunction<CT, FT> createReferenceShapeFunctionRelativeTo(CT cell);
	
	@Override
	ScalarShapeFunctionWithReferenceShapeFunction<CT, FT> createReferenceShapeFunctionRelativeTo(FT face);
	
}
