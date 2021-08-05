package mixed;

import basic.Function;
import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ConcurrentHashMap;

public class MixedFunction implements Function < MixedValue,
	MixedGradient,
	MixedHessian>
{
	private final ScalarFunction pressureFunction;
	private final VectorFunction velocityFunction;
	private boolean overridesValue;
	
	public MixedFunction()
	{
		overridesValue = true;
		pressureFunction = null;
		velocityFunction = null;
	}
	public MixedFunction(ScalarFunction pressureFunction, VectorFunction velocityFunction)
	{
		overridesValue = false;
		this.pressureFunction = pressureFunction;
		this.velocityFunction = velocityFunction;
	}
	
	public ScalarFunction getPressureFunction()
	{
		if(!hasPressureFunction())
			throw new IllegalStateException("Is not pressure Function");
		return pressureFunction;
	}
	
	public VectorFunction getVelocityFunction()
	{
		if(!hasVelocityFunction())
			throw new IllegalStateException("Is not velocity Function");
		return velocityFunction;
	}
	
	public MixedFunction(@NotNull ScalarFunction pressureFunction)
	{
		this.pressureFunction = pressureFunction;
		this.velocityFunction = null;
	}
	public MixedFunction(@NotNull VectorFunction velocityFunction)
	{
		this.pressureFunction = null;
		this.velocityFunction = velocityFunction;
	}
	public boolean hasPressureFunction()
	{
		return pressureFunction!=null;
	}
	public boolean hasVelocityFunction()
	{
		return velocityFunction!=null;
	}
	@Override
	public int getDomainDimension()
	{
		if(hasPressureFunction())
			return Objects.requireNonNull(getPressureFunction()).getDomainDimension();
		else
			return Objects.requireNonNull(getVelocityFunction()).getDomainDimension();
	}
	
	@Override
	public MixedValue defaultValue()
	{
		if(velocityFunction != null)
			return new MixedValue(getVelocityFunction().getRangeDimension()+1);
		if(pressureFunction != null)
			return new MixedValue(getPressureFunction().getDomainDimension()+1);
		return new MixedValue(0);
	}
	
	@Override
	public MixedGradient defaultGradient()
	{
		if(velocityFunction != null)
			return new MixedGradient(getVelocityFunction().getRangeDimension()+1);
		if(pressureFunction != null)
			return new MixedGradient(getPressureFunction().getDomainDimension()+1);
		return new MixedGradient(0);
	}
	
	@Override
	public MixedHessian defaultHessian()
	{
		if(velocityFunction != null)
			return new MixedHessian(getVelocityFunction().getRangeDimension()+1);
		if(pressureFunction != null)
			return new MixedHessian(getPressureFunction().getDomainDimension()+1);
		return new MixedHessian(0);
	}
	
	@Override
	public MixedValue value(CoordinateVector pos)
	{
		if(hasPressureFunction())
			return new PressureValue(Objects.requireNonNull(getPressureFunction()).value(pos));
		if(hasVelocityFunction())
			return new VelocityValue(Objects.requireNonNull(getVelocityFunction()).value(pos));
		if(overridesValue)
			throw new IllegalStateException("needs to override value");
		throw new IllegalStateException("neither pressure nor velocity function");
	}
	
	@Override
	public MixedGradient gradient(CoordinateVector pos)
	{
		if(hasPressureFunction())
			return new PressureGradient(Objects.requireNonNull(getPressureFunction()).gradient(pos));
		if(hasVelocityFunction())
			return new VelocityGradient(Objects.requireNonNull(getVelocityFunction()).gradient(pos));
		if(overridesValue)
			throw new IllegalStateException("needs to override gradient");
		throw new IllegalStateException("neither pressure nor velocity function");
	}
	public Map<CoordinateVector, Double> pressureValuesInPoints(List<CoordinateVector> points)
	{
		ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point, value(point).getPressure()));
		return ret;
	}
	public Map<CoordinateVector, Double> pressureValuesInPointsAtTime(List<CoordinateVector> points, double t)
	{
		ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point.addTime(t), value(point).getPressure()));
		return ret;
	}
	public Map<CoordinateVector, Double> velocityComponentsInPoints(List<CoordinateVector> points, int component)
	{
		ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point, value(point).getVelocity().at(component)));
		return ret;
	}
	
	public Map<CoordinateVector, Double> velocityComponentsInPointsAtTime(List<CoordinateVector> points, int component,
	                                                                      double t)
	{
		ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point.addTime(t),
			value(point).getVelocity().at(component)));
		return ret;
	}
	public Map<CoordinateVector, CoordinateVector> velocityValuesInPointsAtTime(List<CoordinateVector> points,
	                                                                      double t)
	{
		ConcurrentHashMap<CoordinateVector, CoordinateVector> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point.addTime(t),
			value(point).getVelocity().addTime(t)));
		return ret;
	}
}
