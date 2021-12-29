package mixed;

import basic.Cell;
import basic.Face;
import basic.ScalarShapeFunction;
import basic.VectorShapeFunction;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;

public abstract class ComposeMixedShapeFunction<CT extends Cell<CT, FT>,
	FT extends Face<CT, FT>, PF extends ScalarShapeFunction<CT, FT
	>, VF extends VectorShapeFunction<CT, FT>>
	implements MixedShapeFunction<CT, FT, PF, VF>
{
	
	private final PF pressureFunction;
	private final VF velocityFunction;
	private final MixedNodeFunctional nodeFunctional;
	
	public ComposeMixedShapeFunction(@NotNull final PF pressureFunction)
	{
		this.pressureFunction = pressureFunction;
		this.velocityFunction = null;
		nodeFunctional = MixedNodeFunctional.pressureFunctional(pressureFunction.getNodeFunctional());
	}
	
	public ComposeMixedShapeFunction(@NotNull final VF velocityFunction)
	{
		this.pressureFunction = null;
		this.velocityFunction = velocityFunction;
		nodeFunctional = MixedNodeFunctional.velocityFunctional(velocityFunction.getNodeFunctional());
	}
	
	@Override
	public Collection<CT> getCells()
	{
		if (hasPressureFunction())
			return getPressureFunction().getCells();
		else
			return getVelocityFunction().getCells();
	}
	
	@Override
	public Collection<FT> getFaces()
	{
		if (hasPressureFunction())
			return getPressureFunction().getFaces();
		else
			return getVelocityFunction().getFaces();
	}
	
	@Override
	public boolean hasPressureFunction()
	{
		return pressureFunction != null;
	}
	
	@Override
	public boolean hasVelocityFunction()
	{
		return velocityFunction != null;
	}
	
	@Override
	public abstract boolean equals(Object other);
	
	@Override
	public abstract int hashCode();
	
	@Override
	public PF getPressureFunction()
	{
		return pressureFunction;
	}
	
	@Override
	public VF getVelocityFunction()
	{
		return velocityFunction;
	}
	
	@Override
	public MixedValue valueInCell(final CoordinateVector pos, final CT cell)
	{
		if (hasPressureFunction())
			return new PressureValue(getPressureFunction().valueInCell(pos, cell), pos.getLength());
		else
			return new VelocityValue(getVelocityFunction().valueInCell(pos, cell));
	}
	
	@Override
	public MixedGradient gradientInCell(final CoordinateVector pos, final CT cell)
	{
		if (hasPressureFunction())
			return new PressureGradient(getPressureFunction().gradientInCell(pos, cell));
		else
			return new VelocityGradient(getVelocityFunction().gradientInCell(pos, cell));
	}
	
	@Override
	public MixedNodeFunctional getNodeFunctional()
	{
		return nodeFunctional;
	}
}
