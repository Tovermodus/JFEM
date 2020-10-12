package basic;

import linalg.CoordinateVector;
import linalg.Matrix;

public class IntegralNodeFunctional implements NodeFunctional<ScalarFunction, Double, CoordinateVector, Matrix>
{
	private final CoordinateVector point;
	public IntegralNodeFunctional(CoordinateVector point)
	{
		this.point = point;
	}
	
	public CoordinateVector getPoint()
	{
		return point;
	}
	
	@Override
	public double evaluate(ScalarFunction func)
	{
		 return func.value(point);
	}
}
