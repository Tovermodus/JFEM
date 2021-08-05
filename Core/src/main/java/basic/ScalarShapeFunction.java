package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

import java.util.Map;
import java.util.Set;

public interface ScalarShapeFunction<CT extends Cell<CT,FT>, FT extends Face<CT,FT>> extends ScalarFunction, ShapeFunction<CT
	,FT,
	Double, CoordinateVector, CoordinateMatrix>
{
	@Override
	default <ST extends ShapeFunction<CT, FT, Double, CoordinateVector, CoordinateMatrix>> Map<Integer, Double>
	prolongate(Set<ST> refinedFunctions)
	{
		throw new UnsupportedOperationException("Not Yet Implemented");
	}
	
	@Override
	NodeFunctional<Double, CoordinateVector, CoordinateMatrix> getNodeFunctional();
	
	
	@Override
	default Double value(CoordinateVector pos)
	{
		for(CT cell: getCells())
			if(cell.isInCell(pos))
				return valueInCell(pos, cell);
		return 0.;
	}
	
	@Override
	default CoordinateVector gradient(CoordinateVector pos)
	{
		for(CT cell: getCells())
			if(cell.isInCell(pos))
				return gradientInCell(pos, cell);
		return new CoordinateVector(getDomainDimension());
	}
	@Override
	default CoordinateMatrix hessian(CoordinateVector pos)
	{
		for(CT cell: getCells())
			if(cell.isInCell(pos))
				return hessianInCell(pos, cell);
		return new CoordinateMatrix(getDomainDimension(), getDomainDimension());
	}
	@Override
	default Double jumpInValue(FT face, CoordinateVector pos)
	{
		return valueInCell(pos, face.getNormalUpstreamCell()) - valueInCell(pos,
			face.getNormalDownstreamCell());
	}
	@Override
	default CoordinateVector jumpInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell()).sub(
			gradientInCell(pos, face.getNormalDownstreamCell()));
	}
	
	@Override
	default Double averageInValue(FT face, CoordinateVector pos)
	{
		return  0.5*(valueInCell(pos,face.getNormalUpstreamCell())+valueInCell(pos,
			face.getNormalDownstreamCell()));
	}
	
	@Override
	default CoordinateVector averageInDerivative(FT face, CoordinateVector pos)
	{
		return gradientInCell(pos,face.getNormalUpstreamCell()).add(
			gradientInCell(pos,face.getNormalDownstreamCell())).mul(0.5);
	}
	
	@Override
	default CoordinateVector normalAverageInValue(FT face, CoordinateVector pos)
	{
		return face.getNormal().value(pos).mul(0.5*jumpInValue(face,pos));
	}
	
	@Override
	default Double normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		return face.getNormal().value(pos).inner(jumpInDerivative(face, pos));
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
