package basic;

import linalg.*;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

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
	public Map<CoordinateVector, CoordinateVector> valuesInPoints(List<CoordinateVector> points)
	{
		ConcurrentHashMap<CoordinateVector, CoordinateVector> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point, value(point)));
		return ret;
	}
	public Map<CoordinateVector, Double> componentValuesInPoints(List<CoordinateVector> points,
	                                                                       int component)
	{
		ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point, value(point).at(component)));
		return ret;
	}
	public double divergence(CoordinateVector pos)
	{
		double ret = 0;
		for(int i = 0; i < getDomainDimension(); i++)
			ret += gradient(pos).at(i,i);
		return ret;
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
