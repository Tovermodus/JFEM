package basic;

import linalg.*;

public abstract class VectorFunction implements Function<CoordinateVector, CoordinateMatrix, Tensor>
{
	private double getComponentValue(CoordinateVector pos, int component)
	{
		return value(pos).at(component);
	}
	private CoordinateVector getComponentGradient(CoordinateVector pos, int component)
	{
		return (CoordinateVector) gradient(pos).unfoldDimension(0).get(component);
	}
	public ScalarFunction getComponentFunction(int component)
	{
		int domainDimension = getDomainDimension();
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return domainDimension;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return getComponentValue(pos,component);
			}
			
			@Override
			public CoordinateVector gradient(CoordinateVector pos)
			{
				return getComponentGradient(pos, component);
			}
		};
	}
}
