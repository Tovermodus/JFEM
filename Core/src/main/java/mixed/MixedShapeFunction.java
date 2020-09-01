package mixed;

import basic.*;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Map;
import java.util.Set;

public abstract class MixedShapeFunction<CT extends Cell<CT,FT>,
	FT extends Face<CT,FT>, ST extends MixedShapeFunction<CT,FT,ST, PF, VF>, PF extends ScalarShapeFunction<CT,FT
	,PF>,
	VF extends VectorShapeFunction<CT,FT,VF>> extends MixedFunction implements ShapeFunction<CT,
	FT,ST,MixedValue,MixedGradient,MixedHessian>
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
			return getPressureFunction().getCells();
		else
			return getVelocityFunction().getCells();
	}
	@SuppressWarnings("unchecked")
	@Override
	public PF getPressureFunction()
	{
			return (PF) super.getPressureFunction();
	}
	@SuppressWarnings("unchecked")
	@Override
	public VF getVelocityFunction()
	{
		return (VF) super.getVelocityFunction();
	}
	
	@Override
	public Set<FT> getFaces()
	{
		if(isPressure())
			return getPressureFunction().getFaces();
		else
			return getVelocityFunction().getFaces();
	}
	
	@Override
	public NodeFunctional<MixedFunction, MixedValue, MixedGradient, MixedHessian> getNodeFunctional()
	{
		if(isPressure())
			return MixedNodeFunctional.pressureFunctional(getPressureFunction().getNodeFunctional());
		else
			return MixedNodeFunctional.velocityFunctional(getVelocityFunction().getNodeFunctional());
	}
	
	@Override
	public void addFace(FT face)
	{
		if(isPressure())
			getPressureFunction().addFace(face);
		else
			getVelocityFunction().addFace(face);
	}
	
	@Override
	public void addCell(CT cell)
	{
		if(isPressure())
			getPressureFunction().addCell(cell);
		else
			getVelocityFunction().addCell(cell);
	}
	
	@Override
	public MixedValue valueInCell(CoordinateVector pos, CT cell)
	{
		if(isPressure())
			return new PressureValue(getPressureFunction().valueInCell(pos,cell));
		else
			return new VelocityValue(getVelocityFunction().valueInCell(pos,cell));
	}
	
	@Override
	public MixedGradient gradientInCell(CoordinateVector pos, CT cell)
	{
		if(isPressure())
			return new PressureGradient(getPressureFunction().gradientInCell(pos,cell));
		else
			return new VelocityGradient(getVelocityFunction().gradientInCell(pos,cell));
	}
	
	@Override
	public MixedValue jumpInValue(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureValue(getPressureFunction().jumpInValue(face, pos));
		else
			return new VelocityValue(getVelocityFunction().jumpInValue(face, pos));
	}
	
	@Override
	public MixedGradient jumpInDerivative(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureGradient(getPressureFunction().jumpInDerivative(face, pos));
		else
			return new VelocityGradient(getVelocityFunction().jumpInDerivative(face, pos));
	}
	
	@Override
	public MixedValue averageInValue(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureValue(getPressureFunction().averageInValue(face, pos));
		else
			return new VelocityValue(getVelocityFunction().averageInValue(face, pos));
	}
	
	@Override
	public MixedGradient averageInDerivative(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureGradient(getPressureFunction().averageInDerivative(face, pos));
		else
			return new VelocityGradient(getVelocityFunction().averageInDerivative(face, pos));
	}
	
	@Override
	public MixedGradient normalAverageInValue(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureGradient(getPressureFunction().normalAverageInValue(face, pos));
		else
			return new VelocityGradient(getVelocityFunction().normalAverageInValue(face, pos));
	}
	
	@Override
	public MixedValue normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		if(isPressure())
			return new PressureValue(getPressureFunction().normalAverageInDerivative(face, pos));
		else
			return new VelocityValue(getVelocityFunction().normalAverageInDerivative(face, pos));
	}
	
	@Override
	public Map<Integer, Double> prolongate(Set<ST> refinedFunctions)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int compareTo(@NotNull ST o)
	{
		if(o.isPressure() && isPressure())
			return getPressureFunction().compareTo(o.getPressureFunction());
		else if(o.isVelocity() && isVelocity())
			return getVelocityFunction().compareTo(o.getVelocityFunction());
		else if(o.isVelocity() && isPressure())
			return 1;
		else
			return -1;
	}
}
