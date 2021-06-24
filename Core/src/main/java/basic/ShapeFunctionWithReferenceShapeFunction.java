package basic;

import linalg.CoordinateVector;

public interface ShapeFunctionWithReferenceShapeFunction<CT extends CellWithReferenceCell<CT,FT, ET>,
	FT extends FaceWithReferenceFace<CT,FT, ET>,
	ET extends EdgeWithReferenceEdge<CT,FT,ET>, ST extends ShapeFunction<CT,FT,ET, ST,
	valueT,gradientT,hessianT>, valueT,	gradientT,
	hessianT> extends ShapeFunction<CT,FT,ET,ST,valueT,gradientT,hessianT>
{
	public boolean getReferenceShapeFunctionRelativeTo(CellWithReferenceCell<CT,FT,ET> cell,
	                                                         ST otherFunction);
	public boolean getReferenceShapeFunctionRelativeTo(FaceWithReferenceFace<CT,FT,ET> face,
	                                                         ST otherFunction);
	public boolean getReferenceShapeFunctionRelativeTo(EdgeWithReferenceEdge<CT,FT,ET> edge,
	                                                         ST otherFunction);
}
