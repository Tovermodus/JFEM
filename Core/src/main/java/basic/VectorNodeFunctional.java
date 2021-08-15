package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

public class VectorNodeFunctional implements NodeFunctional<CoordinateVector, CoordinateMatrix, CoordinateTensor>
{
	final int component;
	final NodeFunctional<Double, CoordinateVector, CoordinateMatrix> componentNodeFunctional;
	
	public VectorNodeFunctional(final int component, final NodeFunctional<Double, CoordinateVector, CoordinateMatrix> componentNodeFunctional)
	{
		this.component = component;
		this.componentNodeFunctional = componentNodeFunctional;
	}
	
	public double evaluateVectorF(final VectorFunction func)
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
	public double evaluate(final Function<CoordinateVector, CoordinateMatrix, CoordinateTensor> func)
	{
		if (func instanceof VectorFunction)
			evaluateVectorF((VectorFunction) func);
		return evaluateVectorF(VectorFunction.fromRawFunction(func));
	}
	
	@Override
	public boolean usesFace(final Face<?, ?> f)
	{
		return componentNodeFunctional.usesFace(f);
	}
	
	public int getComponent()
	{
		return component;
	}
}
