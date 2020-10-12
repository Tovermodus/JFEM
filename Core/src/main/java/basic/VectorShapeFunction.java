package basic;

import linalg.*;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public abstract class VectorShapeFunction<CT extends Cell<CT,FT, ET>,FT extends Face<CT,FT, ET>, ET extends Edge<CT,
	FT,ET>,
	ST extends VectorShapeFunction<CT,FT,ET,ST>> extends VectorFunction implements ShapeFunction<CT,FT,ET
	,ST,	CoordinateVector,
	CoordinateMatrix, Tensor>, Comparable<ST>
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
	public abstract NodeFunctional<VectorFunction, CoordinateVector,
		CoordinateMatrix, Tensor> getNodeFunctional();
	 
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
	public abstract double divergenceInCell( CoordinateVector pos, CT cell);
	
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
