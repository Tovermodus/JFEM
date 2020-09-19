import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;

import java.util.List;

public class ConvergenceOrderEstimator
{
	public static double estimateL2Vector(List<VectorFunction> functions, List<CoordinateVector> points)
	{
		double average = 0;
		for (int i = 0; i < functions.size()-2; i++)
		{
			double lowerDiff = normL2VecDifference(functions.get(i), functions.get(functions.size()-1), points);
			double higherDiff = normL2VecDifference(functions.get(i+1), functions.get(functions.size()-1),
				points);
			System.out.println(lowerDiff/higherDiff);
			average += lowerDiff/higherDiff;
		}
		return average/(functions.size()-2);
		
	}
	private static double normL2VecDifference(VectorFunction f1, VectorFunction f2, List<CoordinateVector> points)
	{
		double squaredDifferences = 0;
		for(CoordinateVector p: points)
		{
			squaredDifferences += Math.pow(f1.value(p).sub(f2.value(p)).euclidianNorm(),2);
		}
		return Math.sqrt(squaredDifferences);
	}
	public static double estimateL2Scalar(List<ScalarFunction> functions, List<CoordinateVector> points)
	{
		double average = 0;
		for (int i = 0; i < functions.size()-2; i++)
		{
			double lowerDiff = normL2Difference(functions.get(i), functions.get(functions.size()-1), points);
			double higherDiff = normL2Difference(functions.get(i+1), functions.get(functions.size()-1), points);
			System.out.println(lowerDiff/higherDiff);
			average += lowerDiff/higherDiff;
		}
		return average/(functions.size()-2);
		
	}
	private static double normL2Difference(ScalarFunction f1, ScalarFunction f2, List<CoordinateVector> points)
	{
		double squaredDifferences = 0;
		for(CoordinateVector p: points)
		{
			squaredDifferences += Math.pow(f1.value(p)- f2.value(p),2);
		}
		return Math.sqrt(squaredDifferences);
	}
}