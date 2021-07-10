package systems;

import basic.Function;
import basic.FunctionIdentifier;
import basic.PerformanceArguments;
import basic.ScalarFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.netlib.lapack.Dlabad;

import java.util.Map;

public class SystemFunction implements Function<SystemValue,
	SystemGradient,
	SystemHessian>
{
	private final Function<?,?,?>[] functions;
	public SystemFunction(Function<?,?,?>[] functions)
	{
		this.functions = functions;
		if(PerformanceArguments.getInstance().executeChecks)
		{
			for (int i = 0; i < functions.length; i++)
			{
				Function<?, ?, ?> f = functions[i];
				if (SystemParameters.getInstance().signatures[i].getValueT().isInstance(f.defaultValue()))
					throw new IllegalArgumentException("function does not fit signature");
				if (SystemParameters.getInstance().signatures[i].getGradientT().isInstance(f.defaultGradient()))
					throw new IllegalArgumentException("function does not fit signature");
				if (SystemParameters.getInstance().signatures[i].getHessianT().isInstance(f.defaultHessian()))
					throw new IllegalArgumentException("function does not fit signature");
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
		return functions[0].getDomainDimension();
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
		for(int j = 0; j < functions.length; j++)
		for(int i = 0; i < functions.length; i++)
		{
			if(CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[i].getValueT()))
				ret.setComponent((CoordinateVector) functions[i].value(pos), i);
			if(CoordinateMatrix.class.isAssignableFrom(SystemParameters.getInstance().signatures[i].getValueT()))
				ret.setComponent((CoordinateMatrix) functions[i].value(pos), i);
		}
		return ret;
	}
}
