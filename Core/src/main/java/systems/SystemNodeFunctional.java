package systems;

import basic.Function;
import basic.FunctionSignature;
import basic.NodeFunctional;
import basic.PerformanceArguments;

public class SystemNodeFunctional implements NodeFunctional<SystemFunction, SystemValue, SystemGradient, SystemHessian>
{
	private final NodeFunctional<Function<?,?,?>,?,?,?> functional;
	final int component;
	
	public SystemNodeFunctional(NodeFunctional<Function<?,?,?>, ?, ?, ?> functional, int component)
	{
		this.functional = functional;
		this.component = component;
		if (PerformanceArguments.getInstance().executeChecks)
			if(!functional.canEvaluate(SystemParameters.getInstance().signatures[component]))
				throw new IllegalArgumentException("Functional does not fit signatures");
	}
	
	@Override
	public FunctionSignature getFunctionSignature()
	{
		return new FunctionSignature(SystemValue.class, SystemGradient.class, SystemHessian.class);
	}
	
	@Override
	public double evaluate(SystemFunction func)
	{
		return functional.evaluate(func.functions[component]);
	}
}
