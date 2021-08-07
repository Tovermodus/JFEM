package linalg;

import basic.Function;
import basic.VectorFunction;

public class Newton
{
	public static CoordinateVector solve(CoordinateVector initial, CoordinateVector rhs,
	                                                      VectorFunction function)
	{
		CoordinateVector iterate = new CoordinateVector(initial);
		for(int i = 0; i < 100; i++)
		{
			CoordinateVector fx = function.value(iterate).sub(rhs);
			if(fx.almostZero())
			{
				return iterate;
			}
			iterate.subInPlace(function.gradient(iterate).solve(fx));
		}
		System.err.println("Newton did not converge. Final error in maximumnorm was"+ rhs.sub(function.value(iterate)).absMaxElement());
		return iterate;
	}
}
