package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Matrix;
import linalg.Tensor;

public class VectorNodeFunctional implements NodeFunctional<VectorFunction, CoordinateVector, CoordinateMatrix, Tensor>
{
	final int component;
	final NodeFunctional<ScalarFunction, Double, CoordinateVector, Matrix> componentNodeFunctional;
	
	public VectorNodeFunctional(int component, NodeFunctional<ScalarFunction, Double, CoordinateVector, Matrix> componentNodeFunctional)
	{
		this.component = component;
		this.componentNodeFunctional = componentNodeFunctional;
	}
	@Override
	public double evaluate(VectorFunction func)
	{
		return componentNodeFunctional.evaluate(func.getComponentFunction(component));
	}
}
