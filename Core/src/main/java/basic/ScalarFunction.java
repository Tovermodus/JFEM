package basic;

import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Stream;

public interface ScalarFunction extends Function<Double, CoordinateVector, CoordinateMatrix>
{
	static ScalarFunction constantFunction(final double constant)
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return -1;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return constant;
			}
			
			@Override
			public CoordinateVector gradient(final CoordinateVector pos)
			{
				return new CoordinateVector((int) pos.size());
			}
			
			@Override
			public CoordinateMatrix hessian(final CoordinateVector pos)
			{
				return new CoordinateDenseMatrix((int) pos.size(), (int) pos.size());
			}
		};
	}
	
	static ScalarFunction fromRawFunction(final Function<Double, CoordinateVector, CoordinateMatrix> function)
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return function.getDomainDimension();
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return function.value(pos);
			}
			
			@Override
			public CoordinateVector gradient(final CoordinateVector pos)
			{
				return function.gradient(pos);
			}
			
			@Override
			public CoordinateMatrix hessian(final CoordinateVector pos)
			{
				return function.hessian(pos);
			}
		};
	}
	
	default Map<CoordinateVector, Double> valuesInPoints(final List<CoordinateVector> points)
	{
		final ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads) stream = stream.parallel();
		stream.forEach(point -> ret.put(point, value(point)));
		return ret;
	}
	
	default Map<CoordinateVector, Double> valuesInPointsAtTime(final List<CoordinateVector> points, final double t)
	{
		final ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads) stream = stream.parallel();
		stream.forEach(point -> ret.put(point.addCoordinate(t), value(point)));
		return ret;
	}
	
	default VectorFunction getGradientFunction()
	{
		final ScalarFunction me = this;
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return me.getDomainDimension();
			}
			
			@Override
			public int getDomainDimension()
			{
				return me.getDomainDimension();
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return me.gradient(pos);
			}
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return new CoordinateDenseMatrix(me.hessian(pos));
			}
		};
	}
	
	default Double directionalDerivative(final CoordinateVector pos, final CoordinateVector direction)
	{
		return direction.inner(gradient(pos));
	}
	
	default Double delI(final CoordinateVector pos, final int i)
	{
		return gradient(pos).at(i);
	}
	
	default VectorFunction makeIsotropicVectorFunction()
	{
		final ScalarFunction me = this;
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return me.getDomainDimension();
			}
			
			@Override
			public int getDomainDimension()
			{
				return me.getDomainDimension();
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				final CoordinateVector ret = new CoordinateVector(pos.getLength());
				for (int i = 0; i < pos.getLength(); i++)
					ret.set(me.value(pos), i);
				return ret;
			}
		};
	}
	
	@Override
	default Double defaultValue()
	{
		return 0.0;
	}
	
	@Override
	default CoordinateVector defaultGradient()
	{
		return new CoordinateVector(getDomainDimension());
	}
	
	@Override
	default CoordinateMatrix defaultHessian()
	{
		return new CoordinateDenseMatrix(getDomainDimension(), getDomainDimension());
	}
}
