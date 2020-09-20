package basic;

import linalg.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.DoubleStream;

public abstract class ScalarFunction implements Function<Double, CoordinateVector, Matrix>
{
	public static ScalarFunction constantFunction(double constant)
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
	public Map<CoordinateVector, Double> valuesInPoints(List<CoordinateVector> points)
	{
		ConcurrentHashMap<CoordinateVector, Double> ret = new ConcurrentHashMap<>();
		points.stream().parallel().forEach(point->ret.put(point, value(point)));
		return ret;
	}
	public VectorFunction getGradientFunction()
	{
		ScalarFunction me = this;
		return new VectorFunction()
		{
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
	public Double directionalDerivative(CoordinateVector pos, CoordinateVector direction)
	{
		return direction.inner(gradient(pos));
	}
	public Double delI(CoordinateVector pos, int i)
	{
		return gradient(pos).at(i);
	}
	public boolean hasFastEvaluation()
	{
		return false;
	}
	public double fastValue(CoordinateVector pos)
	{
		throw new UnsupportedOperationException();
	}
	public double[] fastGradient(CoordinateVector pos)
	{
		throw new UnsupportedOperationException();
	}
	public double[][] fastHessian(CoordinateVector pos)
	{
		throw new UnsupportedOperationException();
	}
//	default void plot(int pointres, String filename)
//	{
//		double[][] values = new double[pointres][pointres];
//		for (int k = 0; k < pointres; k++)
//		{
//			for (int j = 0; j < pointres; j++)
//			{
//				DoubleTensor pos = DoubleTensor.vectorFromValues(1.0 / (pointres - 1) * k,
//					1.0 / (pointres - 1) * j);
//
//				values[k][j] = value(pos);
//			}
//		}
//		try
//		{
//			BufferedWriter plotWriter = new BufferedWriter(new FileWriter(filename));
//			for (int i = 0; i < pointres; i++)
//			{
//				for (int j = 0; j < pointres; j++)
//				{
//					plotWriter.write(values[i][j] + " ");
//				}
//				plotWriter.newLine();
//			}
//			plotWriter.flush();
//			plotWriter.close();
//		} catch (
//			IOException e)
//		{
//			e.printStackTrace();
//		}
//	}

}
