package systems;

import basic.Function;
import linalg.CoordinateVector;

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
	public SystemValue value(CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
