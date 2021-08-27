package basic;

import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

import java.util.Map;
import java.util.Set;

public interface ScalarShapeFunction<CT extends Cell<CT, FT>, FT extends Face<CT, FT>> extends ScalarFunction, ShapeFunction<CT, FT, Double, CoordinateVector, CoordinateMatrix>
{
	@Override
	default <ST extends ShapeFunction<CT, FT, Double, CoordinateVector, CoordinateMatrix>> Map<Integer, Double> prolongate(final Set<ST> refinedFunctions)
	{
		throw new UnsupportedOperationException("Not Yet Implemented");
	}
	
	@Override
	NodeFunctional<Double, CoordinateVector, CoordinateMatrix> getNodeFunctional();
	
	@Override
	default Double value(final CoordinateVector pos)
	{
		for (final CT cell : getCells())
			if (cell.isInCell(pos)) return valueInCell(pos, cell);
		return 0.;
	}
	
	@Override
	default CoordinateVector gradient(final CoordinateVector pos)
	{
		for (final CT cell : getCells())
			if (cell.isInCell(pos)) return gradientInCell(pos, cell);
		return new CoordinateVector(getDomainDimension());
	}
	
	@Override
	default CoordinateMatrix hessian(final CoordinateVector pos)
	{
		for (final CT cell : getCells())
			if (cell.isInCell(pos)) return hessianInCell(pos, cell);
		return new CoordinateDenseMatrix(getDomainDimension(), getDomainDimension());
	}
	
	@Override
	default Double jumpInValue(final FT face, final CoordinateVector pos)
	{
		return valueInCell(pos, face.getNormalUpstreamCell()) - valueInCell(pos,
		                                                                    face.getNormalDownstreamCell());
	}
	
	@Override
	default CoordinateVector jumpInDerivative(final FT face, final CoordinateVector pos)
	{
		return gradientInCell(pos, face.getNormalUpstreamCell())
			.sub(gradientInCell(pos, face.getNormalDownstreamCell()));
	}
	
	@Override
	default Double averageInValue(final FT face, final CoordinateVector pos)
	{
		return 0.5 * (valueInCell(pos, face.getNormalUpstreamCell()) + valueInCell(pos,
		                                                                           face.getNormalDownstreamCell()));
	}
	
	@Override
	default CoordinateVector averageInDerivative(final FT face, final CoordinateVector pos)
	{
		return gradientInCell(pos, face.getNormalUpstreamCell())
			.add(gradientInCell(pos, face.getNormalDownstreamCell()))
			.mul(0.5);
	}
	
	@Override
	default CoordinateVector normalAverageInValue(final FT face, final CoordinateVector pos)
	{
		return face.getNormal().value(pos).mul(0.5 * jumpInValue(face, pos));
	}
	
	@Override
	default Double normalAverageInDerivative(final FT face, final CoordinateVector pos)
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
