import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;

import java.util.List;

public class ConvergenceOrderEstimator
{
	public static double estimateL2Vector(final List<VectorFunction> functions, final List<CoordinateVector> points)
	{
		double average = 0;
		for (int i = 0; i < functions.size() - 2; i++)
		{
			final double lowerDiff = normL2VecDifference(functions.get(i), functions.get(functions.size() - 1),
			                                             points);
			final double higherDiff = normL2VecDifference(functions.get(i + 1),
			                                              functions.get(functions.size() - 1),
			                                              points);
			System.out.println(lowerDiff + " " + higherDiff);
			System.out.println((1 + i) + ". Estimate of convergence order:" +
				                   Math.log(lowerDiff / higherDiff) / Math.log(2));
			average += Math.log(lowerDiff / higherDiff) / Math.log(2);
		}
		return average / (functions.size() - 2);
	}
	
	public static double normL2VecDifference(final VectorFunction f1, final VectorFunction f2, final List<CoordinateVector> points)
	{
		double squaredDifferences = 0;
		for (final CoordinateVector p : points)
		{
			squaredDifferences += Math.pow(f1.value(p).sub(f2.value(p)).euclidianNorm(), 2);
		}
		return Math.sqrt(squaredDifferences);
	}
	
	public static double estimateL2Scalar(final List<ScalarFunction> functions, final List<CoordinateVector> points)
	{
		double average = 0;
		for (int i = 0; i < functions.size() - 2; i++)
		{
			final double lowerDiff = normL2Difference(functions.get(i), functions.get(functions.size() - 1),
			                                          points);
			final double higherDiff = normL2Difference(functions.get(i + 1), functions.get(functions.size() - 1),
			                                           points);
			System.out.println(lowerDiff + " " + higherDiff);
			System.out.println((1 + i) + "-th Estimate of convergence order:" +
				                   Math.log(lowerDiff / higherDiff) / Math.log(2));
			average += Math.log(lowerDiff / higherDiff) / Math.log(2);
		}
		return average / (functions.size() - 2);
	}
	
	public static double normL2Difference(final ScalarFunction f1, final ScalarFunction f2, final List<CoordinateVector> points)
	{
		double squaredDifferences = 0;
		for (final CoordinateVector p : points)
		{
			squaredDifferences += Math.pow(f1.value(p) - f2.value(p), 2);
		}
		return Math.sqrt(squaredDifferences / points.size());
	}
	
	public static double estimateL20Scalar(final List<ScalarFunction> functions, final List<CoordinateVector> points)
	{
		double average = 0;
		for (int i = 0; i < functions.size() - 2; i++)
		{
			final double lowerDiff = normL20Difference(functions.get(i), functions.get(functions.size() - 1),
			                                           points);
			final double higherDiff = normL20Difference(functions.get(i + 1), functions.get(functions.size() - 1),
			                                            points);
			System.out.println(lowerDiff + " " + higherDiff);
			System.out.println(i + "-th Estimate of convergence order:" +
				                   Math.log(lowerDiff / higherDiff) / Math.log(2));
			average += Math.log(lowerDiff / higherDiff) / Math.log(2);
		}
		return average / (functions.size() - 2);
	}
	
	public static double normL20Difference(final ScalarFunction f1, final ScalarFunction f2, final List<CoordinateVector> points)
	{
		double squaredDifferences = 0;
		double f2Correction = 0;
		for (final CoordinateVector p : points)
			f2Correction += (f1.value(p) - f2.value(p)) / points.size();
		for (final CoordinateVector p : points)
		{
			squaredDifferences += Math.pow(f1.value(p) - (f2Correction + f2.value(p)), 2);
		}
		return Math.sqrt(squaredDifferences);
	}
}
