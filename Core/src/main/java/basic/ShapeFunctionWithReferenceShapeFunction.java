package basic;

import linalg.CoordinateVector;

public interface ShapeFunctionWithReferenceShapeFunction<CT extends CellWithReferenceCell<CT,FT, ET>,
	FT extends FaceWithReferenceFace<CT,FT, ET>,
	ET extends EdgeWithReferenceEdge<CT,FT,ET>, ST extends ShapeFunction<CT,FT,ET, ST,
	valueT,gradientT,hessianT>, valueT,	gradientT,
	hessianT> extends ShapeFunction<CT,FT,ET,ST,valueT,gradientT,hessianT>
{
	
	ST getReferenceShapeFunctionRelativeTo(CT cell);
	ST getReferenceShapeFunctionRelativeTo(FT face);
	default ST getReferenceShapeFunctionRelativeTo(ET edge)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
