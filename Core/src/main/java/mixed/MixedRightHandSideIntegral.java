package mixed;

import basic.*;

public class MixedRightHandSideIntegral<CT extends Cell<CT,FT,ET>,FT extends Face<CT,FT, ET>, ET extends Edge<CT,FT,ET>,
	PF extends ScalarShapeFunction<CT,FT,ET>, VF extends VectorShapeFunction<CT,FT,ET>>
	extends RightHandSideIntegral<CT,MixedShapeFunction<CT,FT,ET,PF,VF>>
{
	private final RightHandSideIntegral<CT, PF> pressureIntegral;
	private final RightHandSideIntegral<CT, VF> velocityIntegral;
	
	private MixedRightHandSideIntegral(RightHandSideIntegral<CT, PF> pressureIntegral,
	                                   RightHandSideIntegral<CT, VF> velocityIntegral)
	{
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	public static <CT extends Cell<CT,FT,ET>,FT extends Face<CT,FT, ET>, ET extends Edge<CT,FT,ET>,
		PF extends ScalarShapeFunction<CT,FT,ET>, VF extends VectorShapeFunction<CT,FT,ET>> MixedRightHandSideIntegral<CT,FT,ET, PF, VF> fromPressureIntegral(RightHandSideIntegral<CT, PF> pressureIntegral)
	{
		return new MixedRightHandSideIntegral<>(pressureIntegral, null);
	}
	
	public static <CT extends Cell<CT,FT,ET>,FT extends Face<CT,FT, ET>, ET extends Edge<CT,FT,ET>,
		PF extends ScalarShapeFunction<CT,FT,ET>, VF extends VectorShapeFunction<CT,FT,ET>> MixedRightHandSideIntegral<CT,FT,ET ,PF,
		VF> fromVelocityIntegral(RightHandSideIntegral<CT
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
	                                            MixedShapeFunction<CT,FT,ET, PF, VF> shapeFunction1)
	{
		if (isPressureIntegral())
		{
			if (shapeFunction1.isVelocity())
				return 0;
			return pressureIntegral.evaluateRightHandSideIntegral(cell, shapeFunction1.getPressureShapeFunction());
		} else
		{
			if (shapeFunction1.isPressure())
				return 0;
			return velocityIntegral.evaluateRightHandSideIntegral(cell, shapeFunction1.getVelocityShapeFunction());
		}
	}
}
