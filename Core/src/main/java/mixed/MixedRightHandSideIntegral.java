package mixed;

import basic.Cell;
import basic.RightHandSideIntegral;
import basic.ScalarShapeFunction;
import basic.VectorShapeFunction;

public class MixedRightHandSideIntegral<CT extends Cell<CT, ?>,
	PF extends ScalarShapeFunction<CT, ?>, VF extends VectorShapeFunction<CT, ?>, MF extends ComposeMixedShapeFunction<CT,
	?, PF, VF>>
	extends RightHandSideIntegral<CT, MF>
{
	private final RightHandSideIntegral<CT, PF> pressureIntegral;
	private final RightHandSideIntegral<CT, VF> velocityIntegral;
	
	private MixedRightHandSideIntegral(final RightHandSideIntegral<CT, PF> pressureIntegral,
	                                   final RightHandSideIntegral<CT, VF> velocityIntegral)
	{
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	public static <CT extends Cell<CT, ?>,
		PF extends ScalarShapeFunction<CT, ?>, VF extends VectorShapeFunction<CT, ?>,
		MF extends ComposeMixedShapeFunction<CT,
			?, PF, VF>> MixedRightHandSideIntegral<CT, PF, VF, MF> fromPressureIntegral(
		final RightHandSideIntegral<CT,
			PF> pressureIntegral)
	{
		return new MixedRightHandSideIntegral<>(pressureIntegral, null);
	}
	
	public static <CT extends Cell<CT, ?>,
		PF extends ScalarShapeFunction<CT, ?>, VF extends VectorShapeFunction<CT, ?>,
		MF extends ComposeMixedShapeFunction<CT,
			?, PF, VF>> MixedRightHandSideIntegral<CT, PF,
		VF, MF> fromVelocityIntegral(final RightHandSideIntegral<CT
		, VF> velocityIntegral)
	{
		return new MixedRightHandSideIntegral<>(null, velocityIntegral);
	}
	
	public boolean isPressureIntegral()
	{
		return this.pressureIntegral != null;
	}
	
	public boolean isVelocityIntegral()
	{
		return this.velocityIntegral != null;
	}
	
	@Override
	public double evaluateRightHandSideIntegral(final CT cell,
	                                            final MF shapeFunction1)
	{
		if (isPressureIntegral())
		{
			if (shapeFunction1.hasVelocityFunction())
				return 0;
			return pressureIntegral.evaluateRightHandSideIntegral(cell,
			                                                      shapeFunction1.getPressureFunction());
		} else
		{
			if (shapeFunction1.hasPressureFunction())
				return 0;
			return velocityIntegral.evaluateRightHandSideIntegral(cell,
			                                                      shapeFunction1.getVelocityFunction());
		}
	}
}
