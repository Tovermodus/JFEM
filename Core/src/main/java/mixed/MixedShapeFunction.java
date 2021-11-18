package mixed;

import basic.*;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

import java.util.Map;
import java.util.Set;

public interface MixedShapeFunction<CT extends Cell<CT, FT>, FT extends Face<CT, FT>, PF extends ScalarShapeFunction<CT, FT
	>, VF extends VectorShapeFunction<CT, FT>>
	extends MixedFunction, ShapeFunction<CT, FT, MixedValue, MixedGradient, MixedHessian>,
	MixedFunctionOnCells<CT, FT>
{
	@Override
	default int getDomainDimension()
	{
		if (hasPressureFunction())
			return getPressureFunction().getDomainDimension();
		else
			return getVelocityFunction().getDomainDimension();
	}
	
	@Override
	PF getPressureFunction();
	
	@Override
	VF getVelocityFunction();
	
	@Override
	default <ST extends ShapeFunction<CT, FT, MixedValue, MixedGradient, MixedHessian>> Map<Integer, Double> prolongate(
		final Set<ST> refinedFunctions)
	{
		throw new UnsupportedOperationException("Not Yet Implemented");
	}
	
	@Override
	default MixedValue value(final CoordinateVector pos)
	{
		for (final CT cell : getCells())
			if (cell.isInCell(pos)) return valueInCell(pos, cell);
		return defaultValue();
	}
	
	@Override
	default MixedGradient gradient(final CoordinateVector pos)
	{
		for (final CT cell : getCells())
			if (cell.isInCell(pos)) return gradientInCell(pos, cell);
		return defaultGradient();
	}
	
	@Override
	default MixedValue jumpInValue(final FT face, final CoordinateVector pos)
	{
		return valueInCell(pos, face.getNormalUpstreamCell())
			.sub(valueInCell(pos, face.getNormalDownstreamCell()));
	}
	
	@Override
	default MixedGradient jumpInDerivative(final FT face, final CoordinateVector pos)
	{
		return gradientInCell(pos, face.getNormalUpstreamCell())
			.sub(gradientInCell(pos, face.getNormalDownstreamCell()));
	}
	
	@Override
	default MixedValue averageInValue(final FT face, final CoordinateVector pos)
	{
		return valueInCell(pos, face.getNormalUpstreamCell())
			.add(valueInCell(pos, face.getNormalDownstreamCell()))
			.mul(0.5);
	}
	
	@Override
	default MixedGradient averageInDerivative(final FT face, final CoordinateVector pos)
	{
		return gradientInCell(pos, face.getNormalUpstreamCell())
			.add(gradientInCell(pos, face.getNormalDownstreamCell()))
			.mul(0.5);
	}
	
	@Override
	NodeFunctional<MixedValue, MixedGradient, MixedHessian> getNodeFunctional();
	
	@Override
	default MixedGradient normalAverageInValue(final FT face, final CoordinateVector pos)
	{
		final MixedValue jump = jumpInValue(face, pos).mul(0.5);
		final CoordinateVector pressureNormalAverageInValue = face.getNormal()
		                                                          .value(pos)
		                                                          .mul(jump.getPressure());
		final CoordinateMatrix velocityNormalAverageInValue = face.getNormal()
		                                                          .value(pos)
		                                                          .outer(jump.getVelocity());
		return new MixedGradient(pressureNormalAverageInValue, velocityNormalAverageInValue);
	}
	
	@Override
	default MixedValue normalAverageInDerivative(final FT face, final CoordinateVector pos)
	{
		final MixedGradient jump = jumpInDerivative(face, pos).mul(0.5);
		final double pressureNormalAverageInDerivative = face.getNormal()
		                                                     .value(pos)
		                                                     .inner(jump.getPressureGradient());
		final CoordinateVector velocityNormalAverageInDerivative = jump.getVelocityGradient()
		                                                               .mvMul(face.getNormal()
		                                                                          .value(pos));
		return new MixedValue(pressureNormalAverageInDerivative, velocityNormalAverageInDerivative);
	}
}
