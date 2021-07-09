package systems;

import basic.Function;
import basic.FunctionIdentifier;
import linalg.CoordinateVector;
import org.netlib.lapack.Dlabad;

import java.util.Map;

public class SystemFunction implements Function< SystemValue,
	SystemGradient,
	SystemHessian>
{
	@Override
	public int getDomainDimension()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public SystemValue defaultValue()
	{
		return new SystemValue();
	}
	
	@Override
	public SystemGradient defaultGradient()
	{
		return new SystemGradient();
	}
	
	@Override
	public SystemHessian defaultHessian()
	{
		return new SystemHessian();
	}
	
	@Override
	public SystemValue value(CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
