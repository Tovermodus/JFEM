package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Tensor;

public class LagrangeVectorNodeFunctional implements NodeFunctional<VectorFunction, CoordinateVector, CoordinateMatrix,
	Tensor>
{
	
	private CoordinateVector point;
	int component;
	public LagrangeVectorNodeFunctional(CoordinateVector point, int component)
	{
		this.point = point;
		this.component = component;
	}
	
	public CoordinateVector getPoint()
	{
		return point;
	}
	public int getComponent()
	{
		return component;
	}
	@Override
	public double evaluate(VectorFunction func)
	{
		return func.value(point).at(component);
	}
}
