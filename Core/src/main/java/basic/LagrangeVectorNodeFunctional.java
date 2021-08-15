package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

public class LagrangeVectorNodeFunctional implements NodeFunctional<CoordinateVector, CoordinateMatrix,
	CoordinateTensor>
{
	
	private final CoordinateVector point;
	int component;
	
	public LagrangeVectorNodeFunctional(final CoordinateVector point, final int component)
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
	public double evaluate(final Function<CoordinateVector, CoordinateMatrix, CoordinateTensor> func)
	{
		return func.value(point).at(component);
	}
	
	@Override
	public boolean usesFace(final Face<?, ?> f)
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
