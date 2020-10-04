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
		CoordinateMatrix grad = gradient(pos);
		for(int i = 0; i < getDomainDimension(); i++)
			ret += grad.at(i,i);
		return ret;
	}
	public CoordinateVector curl(CoordinateVector pos)
	{
		if(getDomainDimension() == 2)
		{
			throw new IllegalStateException("wrong domain dimension");
		}
		CoordinateVector ret = new CoordinateVector(getDomainDimension());
		CoordinateMatrix grad = gradient(pos);
		ret.set(grad.at(2,1)-grad.at(1,2),0);
		ret.set(grad.at(0,2)-grad.at(2,0),1);
		ret.set(grad.at(1,0)-grad.at(0,1),2);
		return ret;
	}
	public ScalarFunction getDivergenceFunction()
	{
		VectorFunction me  = this;
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return me.getDomainDimension();
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return me.divergence(pos);
			}
		};
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
