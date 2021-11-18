package mixed;

import basic.ScalarFunction;
import basic.VectorFunction;

public class ComposedMixedFunction
	implements MixedFunction
{
	private final ScalarFunction pressureFunction;
	private final VectorFunction velocityFunction;
	
	public ComposedMixedFunction(final ScalarFunction pressureFunction, final VectorFunction velocityFunction)
	{
		this.pressureFunction = pressureFunction;
		this.velocityFunction = velocityFunction;
	}
	
	public ComposedMixedFunction(final ScalarFunction pressureFunction)
	{
		this.pressureFunction = pressureFunction;
		this.velocityFunction = ScalarFunction.constantFunction(0)
		                                      .makeIsotropicVectorFunction();
	}
	
	public ComposedMixedFunction(final VectorFunction velocityFunction)
	{
		this.pressureFunction = ScalarFunction.constantFunction(0);
		this.velocityFunction = velocityFunction;
	}
	
	public ComposedMixedFunction()
	{
		this.pressureFunction = ScalarFunction.constantFunction(0);
		this.velocityFunction = ScalarFunction.constantFunction(0)
		                                      .makeIsotropicVectorFunction();
	}
	
	@Override
	public ScalarFunction getPressureFunction()
	{
		return pressureFunction;
	}
	
	@Override
	public VectorFunction getVelocityFunction()
	{
		return velocityFunction;
	}
	
	@Override
	public boolean hasPressureFunction()
	{
		return true;
	}
	
	@Override
	public boolean hasVelocityFunction()
	{
		return true;
	}
}
