package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface ScalarShapeFunctionWithReferenceShapeFunction<CT extends CellWithReferenceCell<CT,FT,ET>,
	FT extends FaceWithReferenceFace<CT,FT,ET>, ET extends EdgeWithReferenceEdge<CT,FT,ET>>
	extends ShapeFunctionWithReferenceShapeFunction<CT,FT,ET,Double, CoordinateVector, CoordinateMatrix>,
		ScalarShapeFunction<CT,FT,ET>
{
	@Override
	ScalarShapeFunctionWithReferenceShapeFunction<CT, FT, ET> getReferenceShapeFunctionRelativeTo(CT cell);
	
	@Override
	ScalarShapeFunctionWithReferenceShapeFunction<CT, FT, ET> getReferenceShapeFunctionRelativeTo(FT face);
	
	@Override
	default ScalarShapeFunctionWithReferenceShapeFunction<CT, FT, ET> getReferenceShapeFunctionRelativeTo(ET edge)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
