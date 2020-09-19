package mixed;

import basic.Function;
import basic.ScalarFunction;
import basic.ScalarShapeFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;

public class MixedFunction implements Function<MixedValue,
	MixedGradient,
	MixedHessian>
{
	private final ScalarFunction pressureFunction;
	private final VectorFunction velocityFunction;
	private boolean overridesValue;
	
	MixedFunction()
	{
		overridesValue = true;
		pressureFunction = null;
		velocityFunction = null;
	}
	
	public ScalarFunction getPressureFunction()
	{
		if(!isPressure())
			throw new IllegalStateException("Is not pressure Function");
		return pressureFunction;
	}
	
	public VectorFunction getVelocityFunction()
	{
		if(!isVelocity())
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
	public boolean isPressure()
	{
		return pressureFunction!=null && !overridesValue;
	}
	public boolean isVelocity()
	{
		return velocityFunction!=null &&!overridesValue;
	}
	@Override
	public int getDomainDimension()
	{
		if(isPressure())
			return pressureFunction.getDomainDimension();
		else
			return velocityFunction.getDomainDimension();
	}
	
	@Override
	public MixedValue value(CoordinateVector pos)
	{
		if(isPressure())
			return new PressureValue(pressureFunction.value(pos));
		if(isVelocity())
			return new VelocityValue(velocityFunction.value(pos));
		if(overridesValue)
			throw new IllegalStateException("needs to override value");
		throw new IllegalStateException("neither pressure nor velocity function");
	}
	
	@Override
	public MixedGradient gradient(CoordinateVector pos)
	{
		if(isPressure())
			return new PressureGradient(pressureFunction.gradient(pos));
		if(isVelocity())
			return new VelocityGradient(velocityFunction.gradient(pos));
		if(overridesValue)
			throw new IllegalStateException("needs to override gradient");
		throw new IllegalStateException("neither pressure nor velocity function");
	}
	public Map<CoordinateVector, Double> pressureValuesInPoints(List<CoordinateVector> points)
	{
		TreeMap<CoordinateVector, Double> ret = new TreeMap<>();
		points.stream().parallel().forEach(point->ret.put(point, value(point).getPressure()));
		return ret;
	}
	public Map<CoordinateVector, Double> velocityComponentsInPoints(List<CoordinateVector> points, int component)
	{
		ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point, value(point).getVelocity().at(component)));
		return ret;
	}
}
