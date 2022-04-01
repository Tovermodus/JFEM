import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;
import mixed.ComposedMixedFunction;
import mixed.MixedFunction;

public class Kovasznay
{
	public static double lamb(final double reyn)
	{
		return reyn / 2 - Math.sqrt(reyn * reyn / 4 + 4 * Math.PI * Math.PI);
	}
	
	public static VectorFunction rightHandSide(final double reyn, final double nu)
	{
		final double l = lamb(reyn);
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				final double elx = Math.exp(l * pos.x());
				final double pi = Math.PI;
				final double e2lx = Math.exp(2 * l * pos.x());
				final double spy = Math.sin(2 * Math.PI * pos.y());
				final double cpy = Math.cos(2 * Math.PI * pos.y());
				final double xComp = -nu * (-l * l * elx * cpy + 4 * elx * pi * pi * cpy)
					- (1 - elx * cpy) * l * elx * cpy
					+ l * elx * elx * spy * spy
					- l * e2lx;
				final double yComp = -nu * (l * l * l * elx * spy / (2 * pi) - 2 * l * elx * pi * spy)
					+ (1 - elx * cpy) * l * l * elx * spy / (2 * pi)
					+ l * l * elx * elx * spy * cpy / (2 * pi);
				return CoordinateVector.fromValues(xComp, yComp);
			}
		};
	}
	
	public static ScalarFunction pressureReferenceSolution(final double reyn)
	{
		final double l = lamb(reyn);
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return -Math.exp(2 * l * pos.x()) / 2;
			}
		};
	}
	
	public static VectorFunction velocityReferenceSolution(final double reyn)
	{
		final double l = lamb(reyn);
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				final double elx = Math.exp(l * pos.x());
				final double spy = Math.sin(2 * Math.PI * pos.y());
				final double cpy = Math.cos(2 * Math.PI * pos.y());
				return CoordinateVector.fromValues(
					1 - elx * cpy,
					l / (2 * Math.PI) * elx * spy);
			}
		};
	}
	
	public static VectorFunction vectorBoundaryValues(final double reyn)
	{
		final VectorFunction referenceSolution = velocityReferenceSolution(reyn);
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				if (pos.x() == -0.5 || pos.x() == 1.5 || pos.y() == 0 || pos.y() == 2)
					return referenceSolution.value(pos);
				return pos.mul(0);
			}
		};
	}
	
	public static MixedFunction mixedReferenceSolution(final double reyn)
	{
		return new ComposedMixedFunction(pressureReferenceSolution(reyn), velocityReferenceSolution(reyn));
	}
}
