package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

public class LagrangeVectorNodeFunctional implements NodeFunctional<CoordinateVector, CoordinateMatrix,
	CoordinateTensor>
{
	
	private CoordinateVector point;
	int component;
	public LagrangeVectorNodeFunctional(CoordinateVector point, int component)
	{
		this.point = point;
		this.component = component;
	}
	
	@Override
	public FunctionSignature getFunctionSignature()
	{
		return new FunctionSignature(CoordinateVector.class, CoordinateMatrix.class, CoordinateTensor.class);
	}
	
	@Override
	public double evaluate(Function<CoordinateVector, CoordinateMatrix, CoordinateTensor> func)
	{
		return func.value(point).at(component);
	}
	
	@Override
	public boolean usesFace(Face<?, ?> f)
	{
		return f.isOnFace(point);
	}
	
	public CoordinateVector getPoint()
	{
		return point;
	}
	public int getComponent()
	{
		return component;
	}
}
