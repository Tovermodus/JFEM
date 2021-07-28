package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

public interface VectorShapeFunctionWithReferenceShapeFunction<CT extends CellWithReferenceCell<CT,FT,ET>,
	FT extends FaceWithReferenceFace<CT,FT,ET>, ET extends EdgeWithReferenceEdge<CT,FT,ET>>
	extends ShapeFunctionWithReferenceShapeFunction<CT,FT,ET,CoordinateVector, CoordinateMatrix, CoordinateTensor>,
		VectorShapeFunction<CT,FT,ET>
{
	@Override
	VectorShapeFunctionWithReferenceShapeFunction<CT, FT, ET> createReferenceShapeFunctionRelativeTo(CT cell);
	
	@Override
	VectorShapeFunctionWithReferenceShapeFunction<CT, FT, ET> createReferenceShapeFunctionRelativeTo(FT face);
	
	@Override
	default VectorShapeFunctionWithReferenceShapeFunction<CT, FT, ET> createReferenceShapeFunctionRelativeTo(ET edge)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
