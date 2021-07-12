package basic;

import linalg.CoordinateVector;

public interface ShapeFunctionWithReferenceShapeFunction<CT extends CellWithReferenceCell<CT,FT, ET>,
	FT extends FaceWithReferenceFace<CT,FT, ET>,
	ET extends EdgeWithReferenceEdge<CT,FT,ET>, valueT,	gradientT,
	hessianT> extends ShapeFunction<CT,FT,ET,valueT,gradientT,hessianT>
{
	
	ShapeFunctionWithReferenceShapeFunction<CT,FT,ET,valueT,gradientT,hessianT> getReferenceShapeFunctionRelativeTo(CT cell);
	ShapeFunctionWithReferenceShapeFunction<CT,FT,ET,valueT,gradientT,hessianT> getReferenceShapeFunctionRelativeTo(FT face);
	default ShapeFunctionWithReferenceShapeFunction<CT,FT,ET,valueT,gradientT,hessianT> getReferenceShapeFunctionRelativeTo(ET edge)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
