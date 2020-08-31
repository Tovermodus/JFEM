package mixed;

import basic.Function;
import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

public class MixedFunction<PF extends ScalarFunction, VF extends VectorFunction> implements Function<MixedValue,
	MixedGradient,
	MixedHessian>
{
	private PF pressureFunction;
	private VF velocityFunction;
	
	public MixedFunction(@NotNull PF pressureFunction)
	{
		this.pressureFunction = pressureFunction;
		this.velocityFunction = null;
	}
	public MixedFunction(@NotNull VF velocityFunction)
	{
		this.pressureFunction = null;
		this.velocityFunction = velocityFunction;
	}
	
	@Override
	public int getDomainDimension()
	{
		return pressureFunction.getDomainDimension();
	}
	
	@Override
	public MixedValue value(CoordinateVector pos)
	{
		if(velocityFunction == null)
			return new PressureValue(pressureFunction.value(pos));
		if(pressureFunction == null)
			return new VelocityValue(velocityFunction.value(pos));
		throw new IllegalStateException("neither pressure nor velocity function");
	}
	
	@Override
	public MixedGradient gradient(CoordinateVector pos)
	{
		if(velocityFunction == null)
			return new PressureGradient(pressureFunction.gradient(pos));
		if(pressureFunction == null)
			return new VelocityGradient(velocityFunction.gradient(pos));
		throw new IllegalStateException("neither pressure nor velocity function");
	}
}
