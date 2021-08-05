package mixed;

import basic.*;
import linalg.*;

public class MixedNodeFunctional implements NodeFunctional<MixedValue, MixedGradient,
	MixedHessian>
{
	private final NodeFunctional<Double, CoordinateVector, CoordinateMatrix> pressureFunctional;
	private final NodeFunctional<CoordinateVector, CoordinateMatrix, CoordinateTensor> velocityFunctional;
	
	private MixedNodeFunctional(NodeFunctional<Double, CoordinateVector, CoordinateMatrix> pressureFunctional,
	                      NodeFunctional<CoordinateVector, CoordinateMatrix, CoordinateTensor> velocityFunctional)
	{
		this.pressureFunctional = pressureFunctional;
		this.velocityFunctional = velocityFunctional;
	}
	
	public  static MixedNodeFunctional pressureFunctional(NodeFunctional<Double,
		CoordinateVector,
		CoordinateMatrix> pressureFunctional)
	{
		return new MixedNodeFunctional(pressureFunctional, null);
	}
	public  static MixedNodeFunctional velocityFunctional(NodeFunctional<CoordinateVector,
		CoordinateMatrix, CoordinateTensor> velocityFunctional)
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
	public FunctionSignature getFunctionSignature()
	{
		return new FunctionSignature(MixedValue.class, MixedGradient.class, MixedHessian.class);
	}
	
	public double evaluateMixedF(MixedFunction func)
	{
		if(func.hasPressureFunction() && isPressureFunctional())
			return pressureFunctional.evaluate(func.getPressureFunction());
		else if(func.hasVelocityFunction() && isVelocityFunctional())
			return velocityFunctional.evaluate(func.getVelocityFunction());
		else
			return 0;
	}
	@Override
	public double evaluate(Function<MixedValue, MixedGradient, MixedHessian> func)
	{
		if(func instanceof MixedFunction)
			return evaluateMixedF((MixedFunction) func);
		throw new IllegalArgumentException("Needs to be performed on MixedFunction");
	}
	
	@Override
	public boolean usesFace(Face<?, ?> f)
	{
		if(isPressureFunctional())
			return pressureFunctional.usesFace(f);
		if(isVelocityFunctional())
			return velocityFunctional.usesFace(f);
		return true;
	}
	
}
