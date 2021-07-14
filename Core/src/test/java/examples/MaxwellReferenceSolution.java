package examples;

import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;
import mixed.MixedFunction;

public class MaxwellReferenceSolution
{
	public static MixedFunction mixedReferenceSolution()
	{
		return new MixedFunction()
		{
			@Override
			public ScalarFunction getPressureFunction()
			{
				return pressureReferenceSolution();
			}
			
			@Override
			public VectorFunction getVelocityFunction()
			{
				return velocityReferenceSolution();
			}
			
			@Override
			public boolean isPressure()
			{
				return true;
			}
			
			@Override
			public boolean isVelocity()
			{
				return true;
			}
		};
	}
	public static VectorFunction rightHandSide()
	{
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 3;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 3;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(0.10e2 * Math.sin(Math.PI * pos.y()) * Math.PI * Math.PI * Math.pow(Math.sin(Math.PI * pos.x()), 0.2e1) * Math.cos(Math.PI * pos.z()) * Math.sin(Math.PI * pos.z()) * Math.cos(Math.PI * pos.y()) - 0.2e1 * Math.sin(Math.PI * pos.y()) * Math.pow(Math.cos(Math.PI * pos.x()), 0.2e1) * Math.PI * Math.PI * Math.cos(Math.PI * pos.z()) * Math.sin(Math.PI * pos.z()) * Math.cos(Math.PI * pos.y())
,
				-0.20e2 * Math.sin(Math.PI * pos.z()) * Math.PI * Math.PI * Math.pow(Math.sin(Math.PI * pos.y()), 0.2e1) * Math.cos(Math.PI * pos.x()) * Math.sin(Math.PI * pos.x()) * Math.cos(Math.PI * pos.z()) + 0.4e1 * Math.sin(Math.PI * pos.z()) * Math.pow(Math.cos(Math.PI * pos.y()), 0.2e1) * Math.PI * Math.PI * Math.cos(Math.PI * pos.x()) * Math.sin(Math.PI * pos.x()) * Math.cos(Math.PI * pos.z())
,
					-0.2e1 * Math.sin(Math.PI * pos.x()) * Math.cos(Math.PI * pos.y()) * Math.sin(Math.PI * pos.y()) * Math.pow(Math.cos(Math.PI * pos.z()), 0.2e1) * Math.PI * Math.PI * Math.cos(Math.PI * pos.x()) + 0.10e2 * Math.pow(Math.sin(Math.PI * pos.z()), 0.2e1) * Math.cos(Math.PI * pos.y()) * Math.sin(Math.PI * pos.y()) * Math.cos(Math.PI * pos.x()) * Math.PI * Math.PI * Math.sin(Math.PI * pos.x())
);
			
			}
		};
	}
	public static ScalarFunction pressureReferenceSolution()
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 3;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return 0.;
			}
			
		};
	}
	public static VectorFunction velocityReferenceSolution()
	{
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 3;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 3;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(
					Math.pow(Math.sin(Math.PI * pos.x()), 0.2e1) * Math.cos(Math.PI * pos.y()) * Math.sin(Math.PI * pos.y()) * Math.cos(Math.PI * pos.z()) * Math.sin(Math.PI * pos.z()),
					-0.2e1 * Math.pow(Math.sin(Math.PI * pos.y()), 0.2e1) * Math.cos(Math.PI * pos.x()) * Math.sin(Math.PI * pos.x()) * Math.cos(Math.PI * pos.z()) * Math.sin(Math.PI * pos.z()),
					Math.pow(Math.sin(Math.PI * pos.z()), 0.2e1) * Math.cos(Math.PI * pos.y()) * Math.sin(Math.PI * pos.y()) * Math.cos(Math.PI * pos.x()) * Math.sin(Math.PI * pos.x()));
			}
			
		};
	}
	public static ScalarFunction pressureBoundaryValues()
	{
		ScalarFunction referenceSolution = pressureReferenceSolution();
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 3;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				if(pos.x() == 0|| pos.x() == 1 || pos.y() == 0 || pos.y() == 1|| pos.z() == 0 || pos.z() == 1)
					return referenceSolution.value(pos);
				return 0.;
			}
		};
	}
	public static VectorFunction vectorBoundaryValues()
	{
		VectorFunction referenceSolution = velocityReferenceSolution();
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 3;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 3;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				if(pos.x() == 0|| pos.x() == 1 || pos.y() == 0 || pos.y() == 1|| pos.z() == 0 || pos.z() == 1)
					return referenceSolution.value(pos);
				return pos.mul(0);
			}
		};
	}
}
