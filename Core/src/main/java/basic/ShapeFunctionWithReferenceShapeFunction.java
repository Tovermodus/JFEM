package basic;

public interface ShapeFunctionWithReferenceShapeFunction<CT extends CellWithReferenceCell<CT,FT>,
	FT extends FaceWithReferenceFace<CT,FT>,
	valueT,	gradientT,
	hessianT> extends ShapeFunction<CT,FT, valueT,gradientT,hessianT>
{
	ShapeFunctionWithReferenceShapeFunction<CT,FT, valueT,gradientT,hessianT> createReferenceShapeFunctionRelativeTo(CT cell);
	ShapeFunctionWithReferenceShapeFunction<CT,FT, valueT,gradientT,hessianT> createReferenceShapeFunctionRelativeTo(FT face);
}
