package basic;

import linalg.*;

public class VectorNodeFunctional implements NodeFunctional<VectorFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>
{
	final int component;
	final NodeFunctional<ScalarFunction, Double, CoordinateVector, CoordinateMatrix> componentNodeFunctional;
	
	public VectorNodeFunctional(int component, NodeFunctional<ScalarFunction, Double, CoordinateVector, CoordinateMatrix> componentNodeFunctional)
	{
		this.component = component;
		this.componentNodeFunctional = componentNodeFunctional;
	}
	@Override
	public double evaluate(VectorFunction func)
	{
		return componentNodeFunctional.evaluate(func.getComponentFunction(component));
	}
	
	public NodeFunctional<ScalarFunction, Double, CoordinateVector, CoordinateMatrix> getComponentNodeFunctional()
	{
		return componentNodeFunctional;
	}
	
	public int getComponent()
	{
		return component;
	}
}
