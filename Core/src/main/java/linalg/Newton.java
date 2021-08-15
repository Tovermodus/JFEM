package linalg;

import basic.VectorFunction;

public class Newton
{
	public static CoordinateVector solve(final CoordinateVector initial, final CoordinateVector rhs,
	                                     final VectorFunction function, final int maxIterations)
	{
		final CoordinateVector iterate = new CoordinateVector(initial);
		CoordinateVector step;
		for (int i = 0; i < maxIterations; i++)
		{
			final CoordinateVector fx =
				function.value(iterate).sub(rhs);
//			System.out.println(iterate + "\tit");
//			System.out.println(function.gradient(iterate));
//			System.out.println(fx.absMaxElement() + "\tfx");
			if (fx.almostZero())
			{
				return iterate;
			}
			try
			{
				step = function.gradient(iterate).solve(fx);
			} catch (final RuntimeException e)
			{
				System.out.println("Newton moved into nondifferentiable region, stopping");
				return iterate;
			}
			int backtrs = 0;
			while (function.value(iterate.sub(step)).absMaxElement() > fx.absMaxElement() && backtrs < 1)
			{
				//System.out.println(function.value(iterate.sub(step)).absMaxElement() + "\tlarge");
				step.mulInPlace(0.3);
				backtrs++;
			}
			iterate.subInPlace(step);
		}
		//System.err.println("Newton did not converge. Final error in maximumnorm was" +
		//	                   rhs.sub(function.value(iterate)).absMaxElement());
		return iterate;
	}
}
