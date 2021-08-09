package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public interface ScalarFunction extends Function<Double, CoordinateVector, CoordinateMatrix>
{
	static ScalarFunction constantFunction(double constant)
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return -1;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return constant;
			}
			@Override
			public CoordinateVector gradient(CoordinateVector pos)
			{
				return new CoordinateVector((int) pos.size());
			}
			
			@Override
			public CoordinateMatrix hessian(CoordinateVector pos)
			{
				return new CoordinateMatrix((int)pos.size(), (int)pos.size());
			}
			
		};
	}
	static ScalarFunction fromRawFunction(Function<Double, CoordinateVector, CoordinateMatrix> function)
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return function.getDomainDimension();
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return function.value(pos);
			}
			@Override
			public CoordinateVector gradient(CoordinateVector pos)
			{
				return function.gradient(pos);
			}
			
			@Override
			public CoordinateMatrix hessian(CoordinateVector pos)
			{
				return function.hessian(pos);
			}
			
		};
	}
	default Map<CoordinateVector, Double> valuesInPoints(List<CoordinateVector> points)
	{
		ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point, value(point)));
		return ret;
	}
	default Map<CoordinateVector, Double> valuesInPointsAtTime(List<CoordinateVector> points, double t)
	{
		ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point.addCoordinate(t), value(point)));
		return ret;
	}
	default VectorFunction getGradientFunction()
	{
		ScalarFunction me = this;
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
			public CoordinateVector value(CoordinateVector pos)
			{
				return me.gradient(pos);
			}
			@Override
			public CoordinateMatrix gradient(CoordinateVector pos)
			{
				return new CoordinateMatrix(me.hessian(pos));
			}
		};
	}
	default Double directionalDerivative(CoordinateVector pos, CoordinateVector direction)
	{
		return direction.inner(gradient(pos));
	}
	default Double delI(CoordinateVector pos, int i)
	{
		return gradient(pos).at(i);
	}
	
	default VectorFunction makeIsotropicVectorFunction()
	{
		ScalarFunction me = this;
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
			public CoordinateVector value(CoordinateVector pos)
			{
				CoordinateVector ret = new CoordinateVector(pos.getLength());
				for(int i = 0; i < pos.getLength(); i++)
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
		return new CoordinateMatrix(getDomainDimension(), getDomainDimension());
	}
}
