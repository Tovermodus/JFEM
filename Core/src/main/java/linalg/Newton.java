package linalg;

import basic.VectorFunction;

public class Newton
{
	/*
	f(x_0) + (x-x_0)^T Grad f(x_0) = b
	(x-x_0)^T Grad f(x_0) = b-f(x_0)
	(x-x_0) = (b - f(x_0))^T Grad f(x_0)^{-1}
	x = x_0 + (b - f(x_0))^T Grad f(x_0)^{-1}
	 */
	
	public static CoordinateVector solve(final CoordinateVector initial,
	                                     final CoordinateVector rhs,
	                                     final VectorFunction function,
	                                     final int maxIterations)
	{
		final CoordinateVector iterate = new CoordinateVector(initial);
		CoordinateVector step;
		for (int i = 0; i < maxIterations; i++)
		{
			final CoordinateVector bMinFx = rhs.sub(function.value(iterate));
			if (bMinFx.almostZero())
			{
				return iterate;
			}
			try
			{
				final CoordinateMatrix fgrad = function.gradient(iterate);
				final CoordinateDenseMatrix fGradInverse;
				if (fgrad instanceof CoordinateDenseMatrix)
					fGradInverse = ((CoordinateDenseMatrix) fgrad).inverse();
				else fGradInverse = (new CoordinateDenseMatrix(fgrad)).inverse();
				step = fGradInverse.tvMul(bMinFx);
			} catch (final RuntimeException e)
			{
				System.out.println("Newton moved into nondifferentiable region, stopping");
				return iterate;
			}
			int backtrs = 0;
			while (function.value(iterate.add(step))
			               .absMaxElement() > bMinFx.absMaxElement() && backtrs < 1)
			{
				step.mulInPlace(0.3);
				backtrs++;
			}
			iterate.addInPlace(step);
		}
		return iterate;
	}
}
