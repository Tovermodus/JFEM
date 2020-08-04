package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Matrix;
import linalg.Vector;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.Map;

public abstract class ScalarShapeFunction<CT extends Cell<CT,FT,ST>, FT extends Face<CT,FT,ST>,
	ST extends ScalarShapeFunction<CT,FT,ST>> extends ScalarFunction implements ShapeFunction<CT
	,FT,ST,
	Double, Vector, Matrix>
{
	protected int globalIndex;
	
	public double fastValueInCell(CoordinateVector pos, CT cell)
	{
		throw new UnsupportedOperationException();
	}
	public double[] fastGradientInCell(CoordinateVector pos, CT cell)
	{
		throw new UnsupportedOperationException();
	}
	public double[][] fastHessianInCell(CoordinateVector pos, CT cell)
	{
		throw new UnsupportedOperationException();
	}
	@Override
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
	public Double jumpInValue(FT face, CoordinateVector pos)
	{
		return valueInCell(pos, face.getNormalUpstreamCell(pos)) - valueInCell(pos,
			face.getNormalDownstreamCell(pos));
	}
	@Override
	public Vector jumpInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell(pos)).sub(
			gradientInCell(pos, face.getNormalDownstreamCell(pos)));
	}
	
	@Override
	public Double averageInValue(FT face, CoordinateVector pos)
	{
		return  0.5*(valueInCell(pos,face.getNormalUpstreamCell(pos))+valueInCell(pos,
			face.getNormalDownstreamCell(pos)));
	}
	
	@Override
	public Vector averageInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell(pos)).add(
			gradientInCell(pos,face.getNormalDownstreamCell(pos))).mul(0.5);
	}
	
	@Override
	public Vector normalAverageInValue(FT face, CoordinateVector pos)
	{
		return face.getNormal().value(pos).mul(0.5*jumpInValue(face,pos));
	}
	
	@Override
	public Double normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		return face.getNormal().value(pos).inner(jumpInDerivative(face, pos));
	}
}
