package basic;

import linalg.*;
import org.jetbrains.annotations.NotNull;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public interface VectorShapeFunction<CT extends Cell<CT,FT, ET>,FT extends Face<CT,FT, ET>, ET extends Edge<CT,
	FT,ET>>
	extends VectorFunction, ShapeFunction<CT,FT,ET
	,	CoordinateVector,
	CoordinateMatrix, CoordinateTensor>
{
	
	@Override
	default <ST extends ShapeFunction<CT, FT, ET, CoordinateVector, CoordinateMatrix, CoordinateTensor>> Map<Integer,
	Double> prolongate(Set<ST> refinedFunctions)
	{
		throw new UnsupportedOperationException("Not Yet Implemented");
	}
	
	@Override
	default CoordinateVector value(CoordinateVector pos)
	{
		for(CT cell: getCells())
			if(cell.isInCell(pos))
				return valueInCell(pos, cell);
		return new CoordinateVector(pos.getLength());
	}
	
	@Override
	default CoordinateMatrix gradient(CoordinateVector pos)
	{
		for(CT cell: getCells())
			if(cell.isInCell(pos))
				return gradientInCell(pos, cell);
		return new CoordinateMatrix(pos.getLength(),pos.getLength());
	}
	
	@Override
	default CoordinateVector jumpInValue(FT face, CoordinateVector pos)
	{
		return valueInCell(pos, face.getNormalUpstreamCell()).sub(valueInCell(pos,
			face.getNormalDownstreamCell()));
	}
	@Override
	default CoordinateMatrix jumpInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell()).sub(
			gradientInCell(pos, face.getNormalDownstreamCell()));
	}
	
	@Override
	default CoordinateVector averageInValue(FT face, CoordinateVector pos)
	{
		return  valueInCell(pos,face.getNormalUpstreamCell()).add(valueInCell(pos,
			face.getNormalDownstreamCell())).mul(0.5);
	}
	
	@Override
	default CoordinateMatrix averageInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell()).add(
			gradientInCell(pos,face.getNormalDownstreamCell())).mul(0.5);
	}
	
	@Override
	NodeFunctional<CoordinateVector,
		CoordinateMatrix, CoordinateTensor> getNodeFunctional();
	 
	@Override
	default CoordinateMatrix normalAverageInValue(FT face, CoordinateVector pos)
	{
		return (CoordinateMatrix) face.getNormal().value(pos).outer(jumpInValue(face,pos).mul(0.5));
	}
	
	@Override
	default CoordinateVector normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		return jumpInDerivative(face, pos).mvMul(face.getNormal().value(pos));
	}
	default double divergenceInCell( CoordinateVector pos, CT cell)
	{
		if(cell.isInCell(pos))
			return divergence(pos);
		else
			return 0;
	}
	
	/*@Override
	default Map<Integer, Double> prolongate(Set<ST> refinedFunctions)
	{
		Map<Integer, Double> ret = new HashMap<>();
		for(ST shapeFunction:refinedFunctions)
		{
			ret.put(shapeFunction.getGlobalIndex(), shapeFunction.getNodeFunctional().evaluate(this));
		}
		return ret;
	}*/
}
