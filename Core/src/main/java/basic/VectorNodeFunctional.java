package basic;

import linalg.*;

public class VectorNodeFunctional implements NodeFunctional<CoordinateVector, CoordinateMatrix, CoordinateTensor>
{
	final int component;
	final NodeFunctional<Double, CoordinateVector, CoordinateMatrix> componentNodeFunctional;
	
	public VectorNodeFunctional(int component, NodeFunctional<Double, CoordinateVector, CoordinateMatrix> componentNodeFunctional)
	{
		this.component = component;
		this.componentNodeFunctional = componentNodeFunctional;
	}
	public double evaluateVectorF(VectorFunction func)
	{
		return componentNodeFunctional.evaluate(func.getComponentFunction(component));
	}
	
	public NodeFunctional<Double, CoordinateVector, CoordinateMatrix> getComponentNodeFunctional()
	{
		return componentNodeFunctional;
	}
	
	@Override
	public FunctionSignature getFunctionSignature()
	{
		return new FunctionSignature(CoordinateVector.class, CoordinateMatrix.class, CoordinateTensor.class);
	}
	
	@Override
	public double evaluate(Function<CoordinateVector, CoordinateMatrix, CoordinateTensor> func)
	{
		if(func instanceof VectorFunction)
			evaluateVectorF((VectorFunction) func);
		return evaluateVectorF(VectorFunction.fromRawFunction(func));
	}
	
	public int getComponent()
	{
		return component;
	}
}
