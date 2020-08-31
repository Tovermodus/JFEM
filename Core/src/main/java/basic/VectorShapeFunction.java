package basic;

import linalg.CoordinateVector;
import linalg.Matrix;
import linalg.Tensor;
import linalg.Vector;

public abstract class VectorShapeFunction<CT extends Cell<CT,FT,ST>,FT extends Face<CT,FT,ST>,
	ST extends VectorShapeFunction<CT,FT,ST>> extends VectorFunction implements ShapeFunction<CT,FT
	,ST,	Vector,
	Matrix, Tensor>
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
	public Vector jumpInValue(FT face, CoordinateVector pos)
	{
		return valueInCell(pos, face.getNormalUpstreamCell(pos)).sub(valueInCell(pos,
			face.getNormalDownstreamCell(pos)));
	}
	@Override
	public Matrix jumpInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell(pos)).sub(
			gradientInCell(pos, face.getNormalDownstreamCell(pos)));
	}
	
	@Override
	public Vector averageInValue(FT face, CoordinateVector pos)
	{
		return  valueInCell(pos,face.getNormalUpstreamCell(pos)).add(valueInCell(pos,
			face.getNormalDownstreamCell(pos))).mul(0.5);
	}
	
	@Override
	public Matrix averageInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell(pos)).add(
			gradientInCell(pos,face.getNormalDownstreamCell(pos))).mul(0.5);
	}
	
	@Override
	public Matrix normalAverageInValue(FT face, CoordinateVector pos)
	{
		return face.getNormal().value(pos).outer(jumpInValue(face,pos).mul(0.5));
	}
	
	@Override
	public Vector normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		return jumpInDerivative(face, pos).mvMul(face.getNormal().value(pos));
	}
	
}
