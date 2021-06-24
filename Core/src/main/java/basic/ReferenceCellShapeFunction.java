package basic;

import linalg.CoordinateVector;

public interface ReferenceCellShapeFunction<CT extends CellWithReferenceCell<CT,FT, ET>,
	FT extends FaceWithReferenceFace<CT,FT, ET>,
	ET extends EdgeWithReferenceEdge<CT,FT,ET>, ST extends ShapeFunction<CT,FT,ET, ST,
	valueT,gradientT,hessianT>, valueT,	gradientT,
	hessianT> extends ShapeFunction<CT,FT,ET,ST,valueT,gradientT,hessianT>
{
	@Override
	valueT value(CoordinateVector pos);
	
	@Override
	default gradientT gradient(CoordinateVector pos)
	{
		return ShapeFunction.super.gradient(pos);
	}
	
	@Override
	default hessianT hessian(CoordinateVector pos)
	{
		return ShapeFunction.super.hessian(pos);
	}
	
	@Override
	valueT valueInCell(CoordinateVector pos, CT cell);
	
	@Override
	default gradientT gradientInCell(CoordinateVector pos, CT cell)
	{
		return ShapeFunction.super.gradientInCell(pos, cell);
	}
	
	@Override
	default hessianT hessianInCell(CoordinateVector pos, CT cell)
	{
		return ShapeFunction.super.hessianInCell(pos, cell);
	}
	
	@Override
	valueT jumpInValue(FT face, CoordinateVector pos);
	
	@Override
	gradientT jumpInDerivative(FT face, CoordinateVector pos);
	
	@Override
	valueT averageInValue(FT face, CoordinateVector pos);
	
	@Override
	gradientT averageInDerivative(FT face, CoordinateVector pos);
	
	@Override
	gradientT normalAverageInValue(FT face, CoordinateVector pos);
	
	@Override
	valueT normalAverageInDerivative(FT face, CoordinateVector pos);
}
