import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.*;

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
			
			@Override
			public CoordinateVector gradient(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(-4*(pos.x()+3)*(pos.y()+1)/((pos.x()+3)*(pos.x()+3)+(pos.y()+1)*(pos.y()+1))/((pos.x()+3)*(pos.x()+3)+(pos.y()+1)*(pos.y()+1)),
					2*(pos.x()*pos.x()+6*pos.x()-pos.y()*pos.y()-2*pos.y()+8)/(pos.x()*pos.x()+6*pos.x()+pos.y()*pos.y()+2*pos.y()+10)/(pos.x()*pos.x()+6*pos.x()+pos.y()*pos.y()+2*pos.y()+10));
			}
			
			@Override
			public DenseMatrix hessian(CoordinateVector pos)
			{
				return DenseMatrix.squareMatrixFromValues(
					16*Math.pow(pos.x()+3,2)*(pos.y()+1)/Math.pow(Math.pow(pos.x()+3,2)+Math.pow(pos.y()+1,2),3)
						- 4*(pos.y()+1)/Math.pow(Math.pow(pos.x()+3,2)+Math.pow(pos.y()+1,2), 2),
					
					16*Math.pow(pos.y()+1,2)*(pos.x()+3)/Math.pow(Math.pow(pos.x()+3,2)+Math.pow(pos.y()+1,2), 3)
						- 4*(pos.x()+3)/Math.pow(Math.pow(pos.x()+3,2)+Math.pow(pos.y()+1,2), 2),
					
					16*Math.pow(pos.y()+1,2)*(pos.x()+3)/Math.pow(Math.pow(pos.x()+3,2)+Math.pow(pos.y()+1,2), 3)
						- 4*(pos.x()+3)/Math.pow(Math.pow(pos.x()+3,2)+Math.pow(pos.y()+1,2), 2),
					
					16*Math.pow(pos.y()+1,3)/Math.pow(Math.pow(pos.x()+3,2)+Math.pow(pos.y()+1,2), 3)
						- 12*(pos.y()+1)/Math.pow(Math.pow(pos.x()+3,2)+Math.pow(pos.y()+1,2), 2));
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
