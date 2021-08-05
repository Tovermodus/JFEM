package mixed;

import basic.*;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public abstract class MixedShapeFunction<CT extends Cell<CT, FT>,
	FT extends Face<CT, FT>, PF extends ScalarShapeFunction<CT, FT
	>, VF extends VectorShapeFunction<CT, FT>>
	extends MixedFunction
	implements ShapeFunction<CT, FT, MixedValue, MixedGradient,
	MixedHessian>, Comparable<MixedShapeFunction<CT,FT, PF,VF>>
{
	
	public MixedShapeFunction(@NotNull PF pressureFunction)
	{
		super(pressureFunction);
	}
	
	public MixedShapeFunction(@NotNull VF velocityFunction)
	{
		super(velocityFunction);
	}
	
	
	@Override
	public Set<CT> getCells()
	{
		if (hasPressureFunction())
			return getPressureShapeFunction().getCells();
		else
			return getVelocityShapeFunction().getCells();
	}
	
	@SuppressWarnings("unchecked")
	
	public PF getPressureShapeFunction()
	{
		return (PF) getPressureFunction();
	}
	
	@SuppressWarnings("unchecked")
	public VF getVelocityShapeFunction()
	{
		return (VF) getVelocityFunction();
	}
	
	@Override
	public Set<FT> getFaces()
	{
		if (hasPressureFunction())
			return getPressureShapeFunction().getFaces();
		else
			return getVelocityShapeFunction().getFaces();
	}
	
	@Override
	public NodeFunctional<MixedValue, MixedGradient, MixedHessian> getNodeFunctional()
	{
		if (hasPressureFunction())
			return MixedNodeFunctional.pressureFunctional(getPressureShapeFunction().getNodeFunctional());
		else
			return MixedNodeFunctional.velocityFunctional(getVelocityShapeFunction().getNodeFunctional());
	}
	
	@Override
	public MixedValue valueInCell(CoordinateVector pos, CT cell)
	{
		if (hasPressureFunction())
			return new PressureValue(getPressureShapeFunction().valueInCell(pos, cell));
		else
			return new VelocityValue(getVelocityShapeFunction().valueInCell(pos, cell));
	}
	
	@Override
	public MixedGradient gradientInCell(CoordinateVector pos, CT cell)
	{
		if (hasPressureFunction())
			return new PressureGradient(getPressureShapeFunction().gradientInCell(pos, cell));
		else
			return new VelocityGradient(getVelocityShapeFunction().gradientInCell(pos, cell));
	}
	
	@Override
	public MixedValue jumpInValue(FT face, CoordinateVector pos)
	{
		if (hasPressureFunction())
			return new PressureValue(getPressureShapeFunction().jumpInValue(face, pos));
		else
			return new VelocityValue(getVelocityShapeFunction().jumpInValue(face, pos));
	}
	
	@Override
	public MixedGradient jumpInDerivative(FT face, CoordinateVector pos)
	{
		if (hasPressureFunction())
			return new PressureGradient(getPressureShapeFunction().jumpInDerivative(face, pos));
		else
			return new VelocityGradient(getVelocityShapeFunction().jumpInDerivative(face, pos));
	}
	
	@Override
	public MixedValue averageInValue(FT face, CoordinateVector pos)
	{
		if (hasPressureFunction())
			return new PressureValue(getPressureShapeFunction().averageInValue(face, pos));
		else
			return new VelocityValue(getVelocityShapeFunction().averageInValue(face, pos));
	}
	
	@Override
	public MixedGradient averageInDerivative(FT face, CoordinateVector pos)
	{
		if (hasPressureFunction())
			return new PressureGradient(getPressureShapeFunction().averageInDerivative(face, pos));
		else
			return new VelocityGradient(getVelocityShapeFunction().averageInDerivative(face, pos));
	}
	
	@Override
	public MixedGradient normalAverageInValue(FT face, CoordinateVector pos)
	{
		if (hasPressureFunction())
			return new PressureGradient(getPressureShapeFunction().normalAverageInValue(face, pos));
		else
			return new VelocityGradient(getVelocityShapeFunction().normalAverageInValue(face, pos));
	}
	
	@Override
	public MixedValue normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		if (hasPressureFunction())
			return new PressureValue(getPressureShapeFunction().normalAverageInDerivative(face, pos));
		else
			return new VelocityValue(getVelocityShapeFunction().normalAverageInDerivative(face, pos));
	}
	
	
	@Override
	public <ST extends ShapeFunction<CT, FT, MixedValue, MixedGradient, MixedHessian>> Map<Integer, Double> prolongate(Set<ST> refinedFunctions)
	{
		Map<Integer, Double> ret = new HashMap<>();
		for (ST shapeFunction : refinedFunctions)
		{
			ret.put(shapeFunction.getGlobalIndex(),
				shapeFunction.getNodeFunctional().evaluate(this));
		}
		return ret;
	}
	
	@SuppressWarnings("unchecked")
	public int compareTo(MixedShapeFunction<CT,FT, PF,VF> o)
	{
		if (o.hasPressureFunction() && hasPressureFunction())
			return ((Comparable<PF>)getPressureShapeFunction()).compareTo(o.getPressureShapeFunction());
		else if (o.hasVelocityFunction() && hasVelocityFunction())
			return ((Comparable<VF>)getVelocityShapeFunction()).compareTo(o.getVelocityShapeFunction());
		else if (o.hasVelocityFunction() && hasPressureFunction())
			return -1;
		else
			return 1;
	}
}
