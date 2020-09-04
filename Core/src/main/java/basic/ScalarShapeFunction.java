package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Matrix;
import linalg.Vector;
import org.jetbrains.annotations.NotNull;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public abstract class ScalarShapeFunction<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	ST extends ScalarShapeFunction<CT,FT,ST>> extends ScalarFunction implements ShapeFunction<CT
	,FT,ST,
	Double, CoordinateVector, Matrix>, Comparable<ST>
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
	public abstract NodeFunctional<ScalarFunction, Double, CoordinateVector, Matrix> getNodeFunctional();
	
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
	public CoordinateVector jumpInDerivative(FT face, CoordinateVector pos)
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
	public CoordinateVector averageInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell(pos)).add(
			gradientInCell(pos,face.getNormalDownstreamCell(pos))).mul(0.5);
	}
	
	@Override
	public CoordinateVector normalAverageInValue(FT face, CoordinateVector pos)
	{
		return face.getNormal().value(pos).mul(0.5*jumpInValue(face,pos));
	}
	
	@Override
	public Double normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		return face.getNormal().value(pos).inner(jumpInDerivative(face, pos));
	}
	
	@Override
	public Map<Integer, Double> prolongate(Set<ST> refinedFunctions)
	{
		Map<Integer, Double> ret = new HashMap<>();
		for(ST shapeFunction:refinedFunctions)
		{
			ret.put(shapeFunction.getGlobalIndex(), shapeFunction.getNodeFunctional().evaluate(this));
		}
		return ret;
	}
}
