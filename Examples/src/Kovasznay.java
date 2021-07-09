import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.DenseMatrix;

public class Kovasznay
{
	static double reynolds = 5;
	static double C = 2*1.70751; // integral over p
	private static double lambda()
	{
		return reynolds/2-Math.sqrt(reynolds*reynolds/4+4*Math.PI*Math.PI);
	}
	//Kovasznay flow
	public static VectorFunction rightHandSide()
	{
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
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(
					1./reynolds*(4*Math.PI*Math.PI-lambda()*lambda())*Math.exp(lambda()*pos.x())*Math.cos(2*Math.PI*pos.y())
					+ lambda()*Math.exp(2*lambda()*pos.x()),
					
					1./reynolds*lambda()/(2*Math.PI)*(4*Math.PI*Math.PI-lambda()*lambda())*Math.exp(lambda()*pos.x())*Math.sin(2*Math.PI*pos.y())
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
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return 0.5*Math.exp(2*lambda()*pos.x())-C;
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
				return 2;
			}
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(1-Math.exp(lambda()* pos.x())*Math.sin(2*Math.PI* pos.x()),
					lambda()/(2*Math.PI)*Math.exp(lambda()* pos.x())*Math.sin(2*Math.PI*pos.y()));
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
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				if(pos.x() == -0.5|| pos.x() == 1.5 || pos.y() == -0.5 || pos.y() == 1.5)
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
				return 2;
			}
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				if(pos.x() == -0.5|| pos.x() == 1.5 || pos.y() == -0.5 || pos.y() == 1.5)
					return referenceSolution.value(pos);
				return pos.mul(0);
			}
		};
	}
}
