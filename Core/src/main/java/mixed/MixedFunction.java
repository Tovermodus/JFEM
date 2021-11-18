package mixed;

import basic.Function;
import basic.PerformanceArguments;
import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Stream;

public interface MixedFunction
	extends Function<MixedValue, MixedGradient, MixedHessian>
{
	
	public ScalarFunction getPressureFunction();
	
	public VectorFunction getVelocityFunction();
	
	@Override
	default int getDomainDimension()
	{
		return getPressureFunction().getDomainDimension();
	}
	
	@Override
	default MixedValue defaultValue()
	{
		return new MixedValue(getDomainDimension());
	}
	
	@Override
	default MixedGradient defaultGradient()
	{
		return new MixedGradient(getDomainDimension());
	}
	
	@Override
	default MixedHessian defaultHessian()
	{
		return new MixedHessian(getDomainDimension());
	}
	
	@Override
	default MixedValue value(final CoordinateVector pos)
	{
		return new MixedValue(getPressureFunction().value(pos), getVelocityFunction().value(pos));
	}
	
	@Override
	default MixedGradient gradient(final CoordinateVector pos)
	{
		final MixedGradient mixedGradient = new MixedGradient(pos.getLength());
		mixedGradient.setPressureGradient(getPressureFunction().gradient(pos));
		mixedGradient.setVelocityGradient(getVelocityFunction().gradient(pos));
		return mixedGradient;
	}
	
	@Override
	default MixedFunction concatenateWith(final VectorFunction f)
	{
		final ScalarFunction concatenatedP = getPressureFunction().concatenateWith(f);
		final VectorFunction concatenatedV = getVelocityFunction().concatenateWith(f);
		return new ComposedMixedFunction(concatenatedP, concatenatedV);
	}
	
	default Map<CoordinateVector, Double> pressureValuesInPoints(final List<CoordinateVector> points)
	{
		final ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(point -> ret.put(point, value(point).getPressure()));
		return ret;
	}
	
	default Map<CoordinateVector, Double> pressureValuesInPointsAtTime(final List<CoordinateVector> points,
	                                                                   final double t)
	{
		final ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(point -> ret.put(point.addCoordinate(t), value(point).getPressure()));
		return ret;
	}
	
	default Map<CoordinateVector, Double> velocityComponentsInPoints(final List<CoordinateVector> points,
	                                                                 final int component)
	{
		final ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(point -> ret.put(point,
		                                value(point).getVelocity()
		                                            .at(component)));
		return ret;
	}
	
	default Map<CoordinateVector, Double> velocityComponentsInPointsAtTime(final List<CoordinateVector> points,
	                                                                       final int component,
	                                                                       final double t)
	{
		final ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(point -> ret.put(point.addCoordinate(t),
		                                value(point).getVelocity()
		                                            .at(component)));
		return ret;
	}
	
	default Map<CoordinateVector, CoordinateVector> velocityValuesInPointsAtTime(final List<CoordinateVector> points,
	                                                                             final double t)
	{
		final ConcurrentHashMap<CoordinateVector, CoordinateVector> ret = new ConcurrentHashMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(point -> ret.put(point.addCoordinate(t),
		                                value(point).getVelocity()
		                                            .addCoordinate(t)));
		return ret;
	}
	
	boolean hasPressureFunction();
	
	boolean hasVelocityFunction();
}
