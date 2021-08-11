package linalg;

import basic.Function;
import basic.VectorFunction;

public class Newton
{
	public static CoordinateVector solve(CoordinateVector initial, CoordinateVector rhs,
	                                                      VectorFunction function)
	{
		CoordinateVector iterate = new CoordinateVector(initial);
		
		for(int i = 0; i < 10; i++)
		{
			CoordinateVector fx = function.value(iterate).sub(rhs);
			if(fx.almostZero())
			{
				return iterate;
			}
			try
			{
				iterate.subInPlace(function.gradient(iterate).solve(fx));
			}
			catch (RuntimeException e)
			{
				//System.out.println("Newton moved into nondifferentiable region, stopping");
				return iterate;
			}
		}
		System.err.println("Newton did not converge. Final error in maximumnorm was"+ rhs.sub(function.value(iterate)).absMaxElement());
		return iterate;
	}
}
