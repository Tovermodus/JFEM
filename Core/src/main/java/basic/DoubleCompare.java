package basic;

public class DoubleCompare
{
	public static boolean almostEqual(final double d1, final double d2)
	{
		return Math.abs(d1 - d2) < PerformanceArguments.getInstance().doubleTolerance * (1 + Math.abs(d1) +
			                                                                                 Math.abs(d2));
	}
	
	public static boolean isInteger(final double d1)
	{
		return Math.round(d1 - Math.round(d1)) <
			PerformanceArguments.getInstance().doubleTolerance * (1 + Math.abs(d1));
	}
	
	public static boolean almostEqual(final double d1, final double d2, final double maxElementInCalculation)
	{
		return Math.abs(d1 - d2) < PerformanceArguments.getInstance().doubleTolerance * (1 + Math.abs(d1) +
			                                                                                 Math.abs(d2) +
			                                                                                 maxElementInCalculation);
	}
	
	public static boolean almostEqualAfterOps(final double d1, final double d2, final long operationNumber)
	{
		return Math.abs(d1 - d2) < PerformanceArguments.getInstance().doubleTolerance * (1 + Math.abs(d1) +
			                                                                                 Math.abs(d2)) *
			operationNumber;
	}
	
	public static boolean almostEqualAfterOps(final double d1, final double d2, final double maxElementInCalculation,
	                                          final long operationNumber)
	{
		return Math.abs(d1 - d2) < PerformanceArguments.getInstance().doubleTolerance * (1 + Math.abs(d1) +
			                                                                                 Math.abs(d2) +
			                                                                                 maxElementInCalculation) *
			operationNumber;
	}
	
	public static boolean almostEqualAfterOps(final double d1, final double d2, final double maxElementInCalculation,
	                                          final long operationNumber, final double extraTol)
	{
		return Math.abs(d1 - d2) <
			(extraTol + PerformanceArguments.getInstance().doubleTolerance) * (1 + Math.abs(d1) + Math.abs(
				d2) + maxElementInCalculation) * operationNumber;
	}
	
	public static int doubleHash(double d)
	{
		d += PerformanceArguments.getInstance().doubleTolerance * 13719 * (Math.abs(d) + 1);
		final double normalized =
			Math.signum(d) * d * d /
				((int) d * d + 1) + 1. / 7 * PerformanceArguments.getInstance().doubleTolerance;
		return (int) (normalized * 2623.7);
	}
}
