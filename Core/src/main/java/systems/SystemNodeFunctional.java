package systems;

import basic.*;

public class SystemNodeFunctional implements NodeFunctional<SystemValue,
	SystemGradient,
	SystemHessian>
{
	private final NodeFunctional functional;
	final int component;
	
	public SystemNodeFunctional(final NodeFunctional<?, ?, ?> functional, final int component)
	{
		this.functional = functional;
		this.component = component;
		if (PerformanceArguments.getInstance().executeChecks)
			if (!functional.canEvaluate(SystemParameters.getInstance().signatures[component]))
				throw new IllegalArgumentException("Functional does not fit signatures");
	}
	
	@Override
	public FunctionSignature getFunctionSignature()
	{
		return new FunctionSignature(SystemValue.class, SystemGradient.class, SystemHessian.class);
	}
	
	@Override
	public double evaluate(final Function<SystemValue, SystemGradient, SystemHessian> func)
	{
		throw new UnsupportedOperationException("Needs SystemFunction");
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public boolean usesFace(final Face<?, ?> f)
	{
		return functional.usesFace(f);
	}
	
	@SuppressWarnings("unchecked")
	public double evaluate(final SystemFunction func)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!functional.canEvaluate(func.getComponentFunction(component).getFunctionSignature()))
				throw new IllegalArgumentException("Functional does not fit signatures");
		return functional.evaluate(func.getComponentFunction(component));
	}
}
