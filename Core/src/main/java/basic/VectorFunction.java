package basic;

import com.google.common.base.Stopwatch;
import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.stream.Stream;

public interface VectorFunction
	extends Function<CoordinateVector, CoordinateMatrix, CoordinateTensor>
{
	private double getComponentValue(final CoordinateVector pos, final int component)
	{
		return value(pos).at(component);
	}
	
	private CoordinateVector getComponentGradient(final CoordinateVector pos, final int component)
	{
		return (CoordinateVector) gradient(pos).unfoldDimension(0)
		                                       .get(component);
	}
	
	int getRangeDimension();
	
	static VectorFunction fromRawFunction(final Function<CoordinateVector, CoordinateMatrix, CoordinateTensor> function)
	{
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return function.defaultValue()
				               .getLength();
			}
			
			@Override
			public int getDomainDimension()
			{
				return function.getDomainDimension();
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return function.value(pos);
			}
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return function.gradient(pos);
			}
			
			@Override
			public CoordinateTensor hessian(final CoordinateVector pos)
			{
				return function.hessian(pos);
			}
		};
	}
	
	@Override
	default CoordinateVector defaultValue()
	{
		return new CoordinateVector(getRangeDimension());
	}
	
	@Override
	default CoordinateMatrix defaultGradient()
	{
		return new CoordinateDenseMatrix(getRangeDimension(), getDomainDimension());
	}
	
	@Override
	default CoordinateTensor defaultHessian()
	{
		return new CoordinateTensor(getRangeDimension(), getDomainDimension(), getDomainDimension());
	}
	
	default Map<CoordinateVector, CoordinateVector> valuesInPoints(final List<CoordinateVector> points)
	{
		final ConcurrentHashMap<CoordinateVector, CoordinateVector> ret = new ConcurrentHashMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads) stream = stream.parallel();
		stream.forEach(point -> ret.put(point, value(point)));
		return ret;
	}
	
	default Map<CoordinateVector, CoordinateVector> valuesInPointsAtTime(final List<CoordinateVector> points,
	                                                                     final double t)
	{
		final ConcurrentSkipListMap<CoordinateVector, CoordinateVector> ret = new ConcurrentSkipListMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads) stream = stream.parallel();
		stream.forEach(point -> ret.put(point.addCoordinate(t), value(point)));
		return ret;
	}
	
	default Map<CoordinateVector, Double> componentValuesInPoints(final List<CoordinateVector> points,
	                                                              final int component)
	{
		final ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads) stream = stream.parallel();
		stream.forEach(point -> ret.put(point, value(point).at(component)));
		return ret;
	}
	
	default double divergence(final CoordinateVector pos)
	{
		double ret = 0;
		final CoordinateMatrix grad = gradient(pos);
		for (int i = 0; i < getDomainDimension(); i++)
			ret += grad.at(i, i);
		return ret;
	}
	
	default CoordinateVector curl(final CoordinateVector pos)
	{
		final Stopwatch s = Stopwatch.createStarted();
		if (getDomainDimension() == 2)
		{
			throw new IllegalStateException("wrong domain dimension");
		}
		final CoordinateVector ret = new CoordinateVector(getDomainDimension());
		final CoordinateMatrix grad = gradient(pos);
		ret.set(grad.at(2, 1) - grad.at(1, 2), 0);
		ret.set(grad.at(0, 2) - grad.at(2, 0), 1);
		ret.set(grad.at(1, 0) - grad.at(0, 1), 2);
		return ret;
	}
	
	default ScalarFunction getDivergenceFunction()
	{
		final VectorFunction me = this;
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return me.getDomainDimension();
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return me.divergence(pos);
			}
		};
	}
	
	default ScalarFunction getComponentFunction(final int component)
	{
		final int domainDimension = getDomainDimension();
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return domainDimension;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return getComponentValue(pos, component);
			}
			
			@Override
			public CoordinateVector gradient(final CoordinateVector pos)
			{
				return getComponentGradient(pos, component);
			}
		};
	}
	
	default CoordinateMatrix jacobian(final CoordinateVector pos)
	{
		return gradient(pos).transpose();
	}
	
	default CoordinateMatrix jacobianTransposed(final CoordinateVector pos)
	{
		return gradient(pos);
	}
	
	@Override
	default VectorFunction concatenateWith(final VectorFunction f)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (f.getRangeDimension() != getDomainDimension())
				throw new IllegalArgumentException("Inner function has the wrong range");
		}
		
		final VectorFunction vectorFunction = this;
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return vectorFunction.getRangeDimension();
			}
			
			@Override
			public int getDomainDimension()
			{
				return vectorFunction.getDomainDimension();
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return vectorFunction.value(f.value(pos));
			}
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return f.gradient(pos)
				        .mmMul(vectorFunction.gradient(f.value(pos)));
			}
		};
	}
}
