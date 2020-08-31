package basic;

import linalg.*;

public abstract class VectorShapeFunction<CT extends Cell<CT,FT>,FT extends Face<CT,FT>,
	ST extends VectorShapeFunction<CT,FT,ST>> extends VectorFunction implements ShapeFunction<CT,FT
	,ST,	CoordinateVector,
	CoordinateMatrix, Tensor>
{
	
	protected int globalIndex;
	
	@
		Override
	public void setGlobalIndex(int index)
	{
		globalIndex = index;
	}
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	@Override
	public CoordinateVector jumpInValue(FT face, CoordinateVector pos)
	{
		return valueInCell(pos, face.getNormalUpstreamCell(pos)).sub(valueInCell(pos,
			face.getNormalDownstreamCell(pos)));
	}
	@Override
	public CoordinateMatrix jumpInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell(pos)).sub(
			gradientInCell(pos, face.getNormalDownstreamCell(pos)));
	}
	
	@Override
	public CoordinateVector averageInValue(FT face, CoordinateVector pos)
	{
		return  valueInCell(pos,face.getNormalUpstreamCell(pos)).add(valueInCell(pos,
			face.getNormalDownstreamCell(pos))).mul(0.5);
	}
	
	@Override
	public CoordinateMatrix averageInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell(pos)).add(
			gradientInCell(pos,face.getNormalDownstreamCell(pos))).mul(0.5);
	}
	
	@Override
	public CoordinateMatrix normalAverageInValue(FT face, CoordinateVector pos)
	{
		return (CoordinateMatrix) face.getNormal().value(pos).outer(jumpInValue(face,pos).mul(0.5));
	}
	
	@Override
	public CoordinateVector normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		return jumpInDerivative(face, pos).mvMul(face.getNormal().value(pos));
	}
	
}
