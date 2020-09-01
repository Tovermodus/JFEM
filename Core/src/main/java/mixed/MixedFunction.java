package mixed;

import basic.Function;
import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

public class MixedFunction implements Function<MixedValue,
	MixedGradient,
	MixedHessian>
{
	private final ScalarFunction pressureFunction;
	private final VectorFunction velocityFunction;
	
	public ScalarFunction getPressureFunction()
	{
		if(isVelocity())
			throw new IllegalStateException("Is not pressure Function");
		return pressureFunction;
	}
	
	public VectorFunction getVelocityFunction()
	{
		if(isPressure())
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
		return pressureFunction!=null;
	}
	public boolean isVelocity()
	{
		return velocityFunction!=null;
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
		throw new IllegalStateException("neither pressure nor velocity function");
	}
	
	@Override
	public MixedGradient gradient(CoordinateVector pos)
	{
		if(isPressure())
			return new PressureGradient(pressureFunction.gradient(pos));
		if(isVelocity())
			return new VelocityGradient(velocityFunction.gradient(pos));
		throw new IllegalStateException("neither pressure nor velocity function");
	}
}
