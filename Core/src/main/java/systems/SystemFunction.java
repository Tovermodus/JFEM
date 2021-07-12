package systems;

import basic.Function;
import basic.FunctionSignature;
import basic.PerformanceArguments;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public class SystemFunction implements Function<SystemValue,
	SystemGradient,
	SystemHessian>
{
	protected final Function<?,?,?>[] functions;
	public final int mainComponent;
	public SystemFunction(Function<?,?,?> function, int component)
	{
		mainComponent = component;
		Function<?,?,?>[] functions = new Function[SystemParameters.getInstance().signatures.length];
		functions[component] = function;
		this.functions = functions;
		if(PerformanceArguments.getInstance().executeChecks)
		{
			if (!SystemParameters.getInstance().signatures[component].getValueT().isInstance(function.defaultValue()))
				throw new IllegalArgumentException("function does not fit signature in value");
			if (!SystemParameters.getInstance().signatures[component].getGradientT().isInstance(function.defaultGradient()))
				throw new IllegalArgumentException("function does not fit signature in " +
					"gradient");
			if (!SystemParameters.getInstance().signatures[component].getHessianT().isInstance(function.defaultHessian()))
				throw new IllegalArgumentException("function does not fit signature in " +
					"hessian");
			if (function.getDomainDimension() != getDomainDimension())
				throw new IllegalArgumentException("Domain Dimensions do not fit");
		}
		
	}
	public SystemFunction(Function<?,?,?>[] functions)
	{
		this.functions = functions;
		mainComponent = 0;
		if(PerformanceArguments.getInstance().executeChecks)
		{
			for (int i = 0; i < functions.length; i++)
			{
				Function<?, ?, ?> f = functions[i];
				if (!SystemParameters.getInstance().signatures[i].getValueT().isInstance(f.defaultValue()))
					throw new IllegalArgumentException("function does not fit signature in value");
				if (!SystemParameters.getInstance().signatures[i].getGradientT().isInstance(f.defaultGradient()))
					throw new IllegalArgumentException("function does not fit signature in " +
						"gradient");
				if (!SystemParameters.getInstance().signatures[i].getHessianT().isInstance(f.defaultHessian()))
					throw new IllegalArgumentException("function does not fit signature in " +
						"hessian");
				if(f.getDomainDimension() != getDomainDimension())
					throw new IllegalArgumentException("Domain Dimensions do not fit");
			}
			if(functions.length == 0)
				throw new IllegalArgumentException("No functions given");
		}
		
	}
	public Function<?,?,?> getComponentFunction(int component)
	{
		return functions[component];
	}
	@Override
	public int getDomainDimension()
	{
		return functions[mainComponent].getDomainDimension();
	}
	
	@Override
	public SystemValue defaultValue()
	{
		return new SystemValue();
	}
	
	@Override
	public SystemGradient defaultGradient()
	{
		return new SystemGradient(getDomainDimension());
	}
	
	@Override
	public SystemHessian defaultHessian()
	{
		return new SystemHessian();
	}
	
	@Override
	public SystemValue value(CoordinateVector pos)
	{
		SystemValue ret = new SystemValue();
		for(int i = 0; i < functions.length; i++)
		{
			if(functions[i] == null)
				continue;
			if(Double.class.isAssignableFrom(SystemParameters.getInstance().signatures[i].getValueT()))
				ret.setComponent((Double) functions[i].value(pos), i);
			if(CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[i].getValueT()))
				ret.setComponent((CoordinateVector) functions[i].value(pos), i);
		}
		return ret;
	}
	@Override
	public SystemGradient gradient(CoordinateVector pos)
	{
		SystemGradient ret = new SystemGradient(getDomainDimension());
		for(int i = 0; i < functions.length; i++)
		{
			if(functions[i] == null)
				continue;
			if(CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[i].getGradientT()))
				ret.setComponent((CoordinateVector) functions[i].gradient(pos), i);
			if(CoordinateMatrix.class.isAssignableFrom(SystemParameters.getInstance().signatures[i].getGradientT()))
				ret.setComponent((CoordinateMatrix) functions[i].gradient(pos), i);
		}
		return ret;
	}
}
