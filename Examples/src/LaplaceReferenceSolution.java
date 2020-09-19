import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;

public class LaplaceReferenceSolution
{
	public static ScalarFunction scalarRightHandSide()
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return 0.;
			}
		};
	}
	public static VectorFunction vectorRightHandSide()
	{
		ScalarFunction scalarRefFunc = scalarRightHandSide();
		return new VectorFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(scalarRefFunc.value(pos),scalarRefFunc.value(pos));
			}
		};
	}
	public static ScalarFunction scalarReferenceSolution()
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return 2*(1+pos.y())/((3+pos.x())*(3+pos.x())+(1+pos.y())*(1+pos.y()));
			}
		};
	}
	public static VectorFunction vectorReferenceSolution()
	{
		ScalarFunction referenceSolution = scalarReferenceSolution();
		return new VectorFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(referenceSolution.value(pos),
					referenceSolution.value(pos));
			}
		};
	}
	public static ScalarFunction scalarBoundaryValues()
	{
		ScalarFunction referenceSolution = scalarReferenceSolution();
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				if(Math.abs(pos.x()) == 1 || Math.abs(pos.y()) == 1)
					return referenceSolution.value(pos);
				return 0.;
			}
		};
	}
	public static VectorFunction vectorBoundaryValues()
	{
		ScalarFunction referenceSolution = scalarBoundaryValues();
		return new VectorFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(referenceSolution.value(pos),
					referenceSolution.value(pos));
			}
		};
	}
}
