package mixed;

import basic.NodeFunctional;
import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.*;

public class MixedNodeFunctional implements NodeFunctional<MixedFunction, MixedValue, MixedGradient,
	MixedHessian>
{
	private final NodeFunctional<ScalarFunction,Double, CoordinateVector, Matrix> pressureFunctional;
	private final NodeFunctional<VectorFunction,CoordinateVector, CoordinateMatrix, Tensor> velocityFunctional;
	
	private MixedNodeFunctional(NodeFunctional<ScalarFunction, Double, CoordinateVector, Matrix> pressureFunctional,
	                      NodeFunctional<VectorFunction, CoordinateVector, CoordinateMatrix, Tensor> velocityFunctional)
	{
		this.pressureFunctional = pressureFunctional;
		this.velocityFunctional = velocityFunctional;
	}
	
	public  static MixedNodeFunctional pressureFunctional(NodeFunctional<ScalarFunction, Double,
		CoordinateVector,
		Matrix> pressureFunctional)
	{
		return new MixedNodeFunctional(pressureFunctional, null);
	}
	public  static MixedNodeFunctional velocityFunctional(NodeFunctional<VectorFunction,
		CoordinateVector,
		CoordinateMatrix, Tensor> velocityFunctional)
	{
		return new MixedNodeFunctional(null,velocityFunctional);
	}
	
	public boolean isPressureFunctional()
	{
		return velocityFunctional == null;
	}
	public boolean isVelocityFunctional()
	{
		return pressureFunctional == null;
	}
	@Override
	public double evaluate(MixedFunction func)
	{
		if(func.isPressure() && isPressureFunctional())
			return pressureFunctional.evaluate(func.getPressureFunction());
		else if(func.isVelocity() && isVelocityFunctional())
			return velocityFunctional.evaluate(func.getVelocityFunction());
		else
			return 0;
	}
}