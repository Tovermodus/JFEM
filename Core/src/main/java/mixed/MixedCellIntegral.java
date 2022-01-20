package mixed;

import basic.*;
import tensorproduct.QuadratureRule1D;

import java.util.Objects;

public class MixedCellIntegral<CT extends Cell<CT, ?>,
	PF extends ScalarShapeFunction<CT, ?>, VF extends VectorShapeFunction<CT, ?>, MF extends ComposeMixedShapeFunction<CT,
	?, PF, VF>>
	extends CellIntegral<CT, MF>
{
	private final CellIntegral<CT, PF> pressureIntegral;
	private final CellIntegral<CT, VF> velocityIntegral;
	
	protected MixedCellIntegral(final CellIntegral<CT, PF> pressureIntegral,
	                            final CellIntegral<CT, VF> velocityIntegral)
	{
		super(null);
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	protected MixedCellIntegral(final Function<?, ?, ?> weight,
	                            final String name,
	                            final QuadratureRule1D quadratureRule1D)
	{
		super(weight, name, quadratureRule1D);
		this.pressureIntegral = null;
		this.velocityIntegral = null;
	}
	
	protected MixedCellIntegral(final String name, final QuadratureRule1D quadratureRule1D)
	{
		super(name, quadratureRule1D);
		this.pressureIntegral = null;
		this.velocityIntegral = null;
	}
	
	protected MixedCellIntegral(final Function<?, ?, ?> weight, final String name)
	{
		this(weight, name, QuadratureRule1D.Gauss5);
	}
	
	protected MixedCellIntegral(final String name)
	{
		this(name, QuadratureRule1D.Gauss5);
	}
	
	public static <CT extends Cell<CT, ?>,
		PF extends ScalarShapeFunction<CT, ?>,
		VF extends VectorShapeFunction<CT, ?>, MF extends ComposeMixedShapeFunction<CT,
		?, PF, VF>> MixedCellIntegral<CT, PF, VF, MF> fromPressureIntegral(final CellIntegral<CT
		, PF> pressureIntegral)
	{
		return new MixedCellIntegral<>(pressureIntegral, null);
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		PF extends ScalarShapeFunction<CT, FT>, VF extends VectorShapeFunction<CT, FT>,
		MF extends ComposeMixedShapeFunction<CT,
			?, PF, VF>> MixedCellIntegral<CT,
		PF,
		VF, MF> fromVelocityIntegral(final CellIntegral<CT
		, VF> velocityIntegral)
	{
		return new MixedCellIntegral<>(null, velocityIntegral);
	}
	
	public boolean isPressureIntegral()
	{
		return pressureIntegral != null;
	}
	
	public boolean isVelocityIntegral()
	{
		return velocityIntegral != null;
	}
	
	public boolean isPressureVelocityIntegral()
	{
		return this.pressureIntegral == null && this.velocityIntegral == null;
	}
	
	protected double evaluatePressureVelocityIntegral(final CT cell,
	                                                  final MF pressureShapeFunction,
	                                                  final MF velocityShapeFunction)
	{
		throw new UnsupportedOperationException("needs to be overwritten");
	}
	
	@Override
	public double evaluateCellIntegral(final CT cell,
	                                   final MF shapeFunction1,
	                                   final MF shapeFunction2)
	{
		if (isPressureIntegral())
		{
			if (shapeFunction1.hasVelocityFunction() || shapeFunction2.hasVelocityFunction())
				return 0;
			return Objects.requireNonNull(pressureIntegral)
			              .evaluateCellIntegral(cell, shapeFunction1.getPressureFunction(),
			                                    shapeFunction2.getPressureFunction());
		} else if (isVelocityIntegral())
		{
			if (shapeFunction1.hasPressureFunction() || shapeFunction2.hasPressureFunction())
				return 0;
			return Objects.requireNonNull(velocityIntegral)
			              .evaluateCellIntegral(cell, shapeFunction1.getVelocityFunction(),
			                                    shapeFunction2.getVelocityFunction());
		} else
		{
			if (shapeFunction1.hasPressureFunction() && shapeFunction2.hasVelocityFunction())
				return evaluatePressureVelocityIntegral(cell, shapeFunction2, shapeFunction1);
			else if (shapeFunction2.hasPressureFunction() && shapeFunction1.hasVelocityFunction())
				return evaluatePressureVelocityIntegral(cell, shapeFunction1, shapeFunction2);
			else
				return 0;
		}
	}
}
