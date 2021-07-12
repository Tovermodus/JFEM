package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Matrix;

public class LagrangeNodeFunctional implements NodeFunctional<Double, CoordinateVector, CoordinateMatrix>
{
	private final CoordinateVector point;
	public LagrangeNodeFunctional(CoordinateVector point)
	{
		this.point = point;
	}
	
	public CoordinateVector getPoint()
	{
		return point;
	}
	
	@Override
	public FunctionSignature getFunctionSignature()
	{
		return new FunctionSignature(Double.class, CoordinateVector.class, CoordinateMatrix.class);
	}
	
	@Override
	public double evaluate(Function<Double, CoordinateVector, CoordinateMatrix> func)
	{
		 return func.value(point);
	}
}
