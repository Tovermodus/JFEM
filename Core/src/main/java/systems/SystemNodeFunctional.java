package systems;

import basic.Function;
import basic.FunctionSignature;
import basic.NodeFunctional;
import basic.PerformanceArguments;
import it.unimi.dsi.fastutil.ints.IntComparator;
import org.jetbrains.annotations.NotNull;
import org.w3c.dom.Node;

public class SystemNodeFunctional implements NodeFunctional<SystemValue,
	SystemGradient,
	SystemHessian>
{
	private final NodeFunctional functional;
	final int component;
	
	public SystemNodeFunctional(NodeFunctional<?,?,?> functional, int component)
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
	public double evaluate(Function<SystemValue, SystemGradient, SystemHessian> func)
	{
		throw new UnsupportedOperationException("Needs SystemFunction");
	}
	
	@SuppressWarnings("unchecked")
	public double evaluate(SystemFunction func)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if(!functional.canEvaluate(func.getComponentFunction(component).getFunctionSignature()))
				throw new IllegalArgumentException("Functional does not fit signatures");
		return functional.evaluate(func.getComponentFunction(component));
	}
	
}
