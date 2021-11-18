package basic;

import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

import java.util.Map;
import java.util.Set;

public interface VectorShapeFunction<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	extends VectorFunctionOnCells<CT, FT>
	, ShapeFunction<CT,
		FT, CoordinateVector, CoordinateMatrix, CoordinateTensor>
{
	
	@Override
	default <ST extends ShapeFunction<CT, FT, CoordinateVector, CoordinateMatrix, CoordinateTensor>> Map<Integer, Double> prolongate(
		final Set<ST> refinedFunctions)
	{
		throw new UnsupportedOperationException("Not Yet Implemented");
	}
	
	@Override
	default CoordinateVector value(final CoordinateVector pos)
	{
		for (final CT cell : getCells())
			if (cell.isInCell(pos)) return valueInCell(pos, cell);
		return new CoordinateVector(pos.getLength());
	}
	
	@Override
	default CoordinateMatrix gradient(final CoordinateVector pos)
	{
		for (final CT cell : getCells())
			if (cell.isInCell(pos)) return gradientInCell(pos, cell);
		return new CoordinateDenseMatrix(pos.getLength(), pos.getLength());
	}
	
	@Override
	default CoordinateVector jumpInValue(final FT face, final CoordinateVector pos)
	{
		return valueInCell(pos, face.getNormalUpstreamCell())
			.sub(valueInCell(pos, face.getNormalDownstreamCell()));
	}
	
	@Override
	default CoordinateMatrix jumpInDerivative(final FT face, final CoordinateVector pos)
	{
		return gradientInCell(pos, face.getNormalUpstreamCell())
			.sub(gradientInCell(pos, face.getNormalDownstreamCell()));
	}
	
	@Override
	default CoordinateVector averageInValue(final FT face, final CoordinateVector pos)
	{
		return valueInCell(pos, face.getNormalUpstreamCell())
			.add(valueInCell(pos, face.getNormalDownstreamCell()))
			.mul(0.5);
	}
	
	@Override
	default CoordinateMatrix averageInDerivative(final FT face, final CoordinateVector pos)
	{
		return gradientInCell(pos, face.getNormalUpstreamCell())
			.add(gradientInCell(pos, face.getNormalDownstreamCell()))
			.mul(0.5);
	}
	
	@Override
	NodeFunctional<CoordinateVector, CoordinateMatrix, CoordinateTensor> getNodeFunctional();
	
	@Override
	default CoordinateMatrix normalAverageInValue(final FT face, final CoordinateVector pos)
	{
		return (CoordinateMatrix) face.getNormal()
		                              .value(pos)
		                              .outer(jumpInValue(face, pos).mul(0.5));
	}
	
	@Override
	default CoordinateVector normalAverageInDerivative(final FT face, final CoordinateVector pos)
	{
		return jumpInDerivative(face, pos).mvMul(face.getNormal()
		                                             .value(pos));
	}
	
	default double divergenceInCell(final CoordinateVector pos, final CT cell)
	{
		if (cell.isInCell(pos)) return divergence(pos);
		else return 0;
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
