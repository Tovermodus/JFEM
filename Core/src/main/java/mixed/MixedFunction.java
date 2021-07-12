package mixed;

import basic.Function;
import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateMatrix;
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
	
	MixedFunction()
	{
		overridesValue = true;
		pressureFunction = null;
		velocityFunction = null;
	}
	
	static MixedFunction fromRawFunction(Function<MixedValue, MixedGradient, MixedHessian> function)
	{
		return new MixedFunction(){
			@Override
			public MixedValue value(CoordinateVector pos)
			{
				return function.value(pos);
			}
			
			@Override
			public MixedGradient gradient(CoordinateVector pos)
			{
				return function.gradient(pos);
			}
			
			@Override
			public MixedValue defaultValue()
			{
				return function.defaultValue();
			}
			
			@Override
			public MixedGradient defaultGradient()
			{
				return function.defaultGradient();
			}
			
			@Override
			public int getDomainDimension()
			{
				return function.getDomainDimension();
			}
		};
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
			return Objects.requireNonNull(pressureFunction).getDomainDimension();
		else
			return Objects.requireNonNull(velocityFunction).getDomainDimension();
	}
	
	@Override
	public MixedValue defaultValue()
	{
		if(velocityFunction != null)
			return new MixedValue(velocityFunction.getRangeDimension()+1);
		if(pressureFunction != null)
			return new MixedValue(pressureFunction.getDomainDimension()+1);
		return new MixedValue(0);
	}
	
	@Override
	public MixedGradient defaultGradient()
	{
		if(velocityFunction != null)
			return new MixedGradient(velocityFunction.getRangeDimension()+1);
		if(pressureFunction != null)
			return new MixedGradient(pressureFunction.getDomainDimension()+1);
		return new MixedGradient(0);
	}
	
	@Override
	public MixedHessian defaultHessian()
	{
		if(velocityFunction != null)
			return new MixedHessian(velocityFunction.getRangeDimension()+1);
		if(pressureFunction != null)
			return new MixedHessian(pressureFunction.getDomainDimension()+1);
		return new MixedHessian(0);
	}
	
	@Override
	public MixedValue value(CoordinateVector pos)
	{
		if(isPressure())
			return new PressureValue(Objects.requireNonNull(pressureFunction).value(pos));
		if(isVelocity())
			return new VelocityValue(Objects.requireNonNull(velocityFunction).value(pos));
		if(overridesValue)
			throw new IllegalStateException("needs to override value");
		throw new IllegalStateException("neither pressure nor velocity function");
	}
	
	@Override
	public MixedGradient gradient(CoordinateVector pos)
	{
		if(isPressure())
			return new PressureGradient(Objects.requireNonNull(pressureFunction).gradient(pos));
		if(isVelocity())
			return new VelocityGradient(Objects.requireNonNull(velocityFunction).gradient(pos));
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
