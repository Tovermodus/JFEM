import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;
import static java.lang.Math.*;

public class StokesReferenceSolution
{
	static double reynolds = 1;
	static double C = 2*1.70751; // integral over p
	public static VectorFunction rightHandSide()
	{
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
				return CoordinateVector.fromValues(
					2*PI*PI*cos(PI*pos.x())*cos(PI*pos.x())*sin(PI*pos.y())*cos(PI*pos.y())
						-6*PI*PI*sin(PI*pos.x())*sin(PI*pos.x())*sin(PI*pos.y())*cos(PI*pos.y())+Math.PI * Math.cos(Math.PI * pos.x()) * Math.sin(Math.PI * pos.y()),
					-2*PI*PI*cos(PI*pos.y())*cos(PI*pos.y())*sin(PI*pos.x())*cos(PI*pos.x())
						+6*PI*PI*sin(PI*pos.y())*sin(PI*pos.y())*sin(PI*pos.x())*cos(PI*pos.x())+Math.sin(Math.PI * pos.x()) * Math.PI * Math.cos(Math.PI * pos.y()));
				
//					16 * (-1./reynolds * pow(PI,2) *  (pow(cos(PI * pos.x()) , 2)-3./4)*sin(PI*pos.y())+sin(PI*pos.x())/16)*PI*cos(PI*pos.y()),
//					16 * (-1./reynolds * pow(PI,2) *  (pow(cos(PI * pos.y()) , 2)-3./4)*sin
//					(PI*pos.x())+sin(PI*pos.y())/16)*PI*cos(PI*pos.x()));
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
				return Math.cos(Math.PI* pos.x())*Math.cos(Math.PI*pos.y());
			}
			
		};
	}
	public static VectorFunction velocityReferenceSolution()
	{
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
				return CoordinateVector.fromValues(
					2*Math.PI*Math.sin(Math.PI*pos.x())*Math.sin(Math.PI*pos.x())*Math.cos(Math.PI*pos.y())*Math.sin(Math.PI*pos.y()),
						2*Math.PI*Math.sin(Math.PI*pos.x())*Math.sin(Math.PI*pos.y())*Math.cos(Math.PI*pos.x())*Math.sin(Math.PI*pos.y())
						);
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
				if(pos.x() == 0|| pos.x() == 1 || pos.y() == 0 || pos.y() == 1)
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
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				if(pos.x() == 0|| pos.x() == 1 || pos.y() == 0 || pos.y() == 1)
					return referenceSolution.value(pos);
				return pos.mul(0);
			}
		};
	}
}
