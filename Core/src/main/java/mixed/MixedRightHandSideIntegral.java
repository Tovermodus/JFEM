package mixed;

import basic.*;

public class MixedRightHandSideIntegral<CT extends Cell<CT,?>,
	PF extends ScalarShapeFunction<CT,?>, VF extends VectorShapeFunction<CT,?>, MF extends MixedShapeFunction<CT,
	?,PF,VF>>
	extends RightHandSideIntegral<CT,MF>
{
	private final RightHandSideIntegral<CT, PF> pressureIntegral;
	private final RightHandSideIntegral<CT, VF> velocityIntegral;
	
	private MixedRightHandSideIntegral(RightHandSideIntegral<CT, PF> pressureIntegral,
	                                   RightHandSideIntegral<CT, VF> velocityIntegral)
	{
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	public static <CT extends Cell<CT, ?>,
		PF extends ScalarShapeFunction<CT, ?>, VF extends VectorShapeFunction<CT, ?>,
		MF extends MixedShapeFunction<CT,
			?, PF, VF>> MixedRightHandSideIntegral<CT, PF, VF, MF> fromPressureIntegral(RightHandSideIntegral<CT,
		PF> pressureIntegral)
	{
		return new MixedRightHandSideIntegral<>(pressureIntegral, null);
	}
	
	public static <CT extends Cell<CT, ?>,
		PF extends ScalarShapeFunction<CT, ?>, VF extends VectorShapeFunction<CT, ?>,
		MF extends MixedShapeFunction<CT,
			?, PF, VF>> MixedRightHandSideIntegral<CT, PF,
		VF, MF> fromVelocityIntegral(RightHandSideIntegral<CT
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
	
	public double evaluateRightHandSideIntegral(CT cell,
	                                            MF shapeFunction1)
	{
		if (isPressureIntegral())
		{
			if (shapeFunction1.hasVelocityFunction())
				return 0;
			return pressureIntegral.evaluateRightHandSideIntegral(cell, shapeFunction1.getPressureShapeFunction());
		} else
		{
			if (shapeFunction1.hasPressureFunction())
				return 0;
			return velocityIntegral.evaluateRightHandSideIntegral(cell, shapeFunction1.getVelocityShapeFunction());
		}
	}
}
