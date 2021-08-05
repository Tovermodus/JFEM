package basic;

public class DoubleCompare
{
	public static boolean almostEqual(double d1, double d2)
	{
		return Math.abs(d1 - d2) < PerformanceArguments.getInstance().doubleTolerance*(1 + Math.abs(d1) + Math.abs(d2));
	}
	public static boolean isInteger(double d1)
	{
		return Math.round(d1 - Math.round(d1)) < PerformanceArguments.getInstance().doubleTolerance*(1+Math.abs(d1));
	}
	public static boolean almostEqual(double d1, double d2, double maxElementInCalculation)
	{
		return Math.abs(d1 - d2) < PerformanceArguments.getInstance().doubleTolerance*(1 + Math.abs(d1) + Math.abs(d2) + maxElementInCalculation);
	}
	public static boolean almostEqualAfterOps(double d1, double d2, long operationNumber)
	{
		return Math.abs(d1 - d2) < PerformanceArguments.getInstance().doubleTolerance*(1 + Math.abs(d1) + Math.abs(d2)) * operationNumber;
	}
	public static boolean almostEqualAfterOps(double d1, double d2, double maxElementInCalculation,
	                                          long operationNumber)
	{
		return Math.abs(d1 - d2) < PerformanceArguments.getInstance().doubleTolerance*(1 + Math.abs(d1) + Math.abs(d2) + maxElementInCalculation) * operationNumber;
	}
	public static boolean almostEqualAfterOps(double d1, double d2, double maxElementInCalculation,
	                                          long operationNumber, double extraTol)
	{
		return Math.abs(d1 - d2) < (extraTol + PerformanceArguments.getInstance().doubleTolerance)*(1 + Math.abs(d1) + Math.abs(d2) + maxElementInCalculation) * operationNumber;
	}
	public static int doubleHash(double d)
	{
		double normalized = d*d/((int)d*d+1);
		return (int)(normalized/(1e4*PerformanceArguments.getInstance().doubleTolerance));
	}
}
