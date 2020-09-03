package mixed;

import basic.*;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Map;
import java.util.Set;

public abstract class MixedShapeFunction<CT extends Cell<CT,FT>,
	FT extends Face<CT,FT>, PF extends ScalarShapeFunction<CT,FT
	,PF> ,
	VF extends VectorShapeFunction<CT,FT,VF>> extends MixedFunction implements ShapeFunction<CT,
	FT,MixedShapeFunction<CT,FT,PF,VF>,MixedValue,MixedGradient,MixedHessian>, Comparable<MixedShapeFunction<CT,
	FT,PF,VF>>
{
	protected int globalIndex;
	
	public MixedShapeFunction(@NotNull PF pressureFunction)
	{
		super(pressureFunction);
	}
	
	public MixedShapeFunction(@NotNull VF velocityFunction)
	{
		super(velocityFunction);
	}
	
	
	@Override
	public void setGlobalIndex(int index)
	{
		globalIndex = index;
	}
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	@Override
	public Set<CT> getCells()
	{
		if(isPressure())
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
		if(isPressure())
			return getPressureShapeFunction().getFaces();
		else
			return getVelocityShapeFunction().getFaces();
	}
	
	@Override
	public NodeFunctional<MixedFunction, MixedValue, MixedGradient, MixedHessian> getNodeFunctional()
	{
		if(isPressure())
			return MixedNodeFunctional.pressureFunctional(getPressureShapeFunction().getNodeFunctional());
		else
			return MixedNodeFunctional.velocityFunctional(getVelocityShapeFunction().getNodeFunctional());
	}
	
	@Override
	public void addFace(FT face)
	{
		if(isPressure())
			getPressureShapeFunction().addFace(face);
		else
			getVelocityShapeFunction().addFace(face);
	}
	
	@Override
	public void addCell(CT cell)
	{
		if(isPressure())
			getPressureShapeFunction().addCell(cell);
		else
			getVelocityShapeFunction().addCell(cell);
	}
	
	@Override
	public MixedValue valueInCell(CoordinateVector pos, CT cell)
	{
		if(isPressure())
			return new PressureValue(getPressureShapeFunction().valueInCell(pos,cell));
		else
			return new VelocityValue(getVelocityShapeFunction().valueInCell(pos,cell));
	}
	
	@Override
	public MixedGradient gradientInCell(CoordinateVector pos, CT cell)
	{
		if(isPressure())
			return new PressureGradient(getPressureShapeFunction().gradientInCell(pos,cell));
		else
			return new VelocityGradient(getVelocityShapeFunction().gradientInCell(pos,cell));
	}
	
	@Override
	public MixedValue jumpInValue(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureValue(getPressureShapeFunction().jumpInValue(face, pos));
		else
			return new VelocityValue(getVelocityShapeFunction().jumpInValue(face, pos));
	}
	
	@Override
	public MixedGradient jumpInDerivative(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureGradient(getPressureShapeFunction().jumpInDerivative(face, pos));
		else
			return new VelocityGradient(getVelocityShapeFunction().jumpInDerivative(face, pos));
	}
	
	@Override
	public MixedValue averageInValue(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureValue(getPressureShapeFunction().averageInValue(face, pos));
		else
			return new VelocityValue(getVelocityShapeFunction().averageInValue(face, pos));
	}
	
	@Override
	public MixedGradient averageInDerivative(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureGradient(getPressureShapeFunction().averageInDerivative(face, pos));
		else
			return new VelocityGradient(getVelocityShapeFunction().averageInDerivative(face, pos));
	}
	
	@Override
	public MixedGradient normalAverageInValue(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureGradient(getPressureShapeFunction().normalAverageInValue(face, pos));
		else
			return new VelocityGradient(getVelocityShapeFunction().normalAverageInValue(face, pos));
	}
	
	@Override
	public MixedValue normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureValue(getPressureShapeFunction().normalAverageInDerivative(face, pos));
		else
			return new VelocityValue(getVelocityShapeFunction().normalAverageInDerivative(face, pos));
	}
	
	@Override
	public Map<Integer, Double> prolongate(Set<MixedShapeFunction<CT,FT,PF,VF>> refinedFunctions)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int compareTo(@NotNull MixedShapeFunction<CT,FT,PF,VF> o)
	{
		if(o.isPressure() && isPressure())
			return getPressureShapeFunction().compareTo( o.getPressureShapeFunction());
		else if(o.isVelocity() && isVelocity())
			return getVelocityShapeFunction().compareTo( o.getVelocityShapeFunction());
		else if(o.isVelocity() && isPressure())
			return 1;
		else
			return -1;
	}
}
