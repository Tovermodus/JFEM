package mixed;

import basic.*;

public class MixedRightHandSideIntegral<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>> extends RightHandSideIntegral<CT,FT,MixedShapeFunction<CT,FT,PF,VF>>
{
	private final RightHandSideIntegral<CT, FT, PF> pressureIntegral;
	private final RightHandSideIntegral<CT, FT, VF> velocityIntegral;
	
	private MixedRightHandSideIntegral(RightHandSideIntegral<CT, FT, PF> pressureIntegral,
	                                   RightHandSideIntegral<CT, FT, VF> velocityIntegral)
	{
		super();
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>> MixedRightHandSideIntegral<CT, FT, PF, VF> fromPressureIntegral(RightHandSideIntegral<CT, FT
		, PF> pressureIntegral)
	{
		return new MixedRightHandSideIntegral<>(pressureIntegral, null);
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>> MixedRightHandSideIntegral<CT, FT, PF,
		VF> fromVelocityIntegral(RightHandSideIntegral<CT, FT
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
	                                            MixedShapeFunction<CT, FT, PF, VF> shapeFunction1)
	{
		if (isPressureIntegral())
		{
			if(shapeFunction1.isVelocity())
				return 0;
			return pressureIntegral.evaluateRightHandSideIntegral(cell, shapeFunction1.getPressureShapeFunction());
		}
		else
		{
			if(shapeFunction1.isPressure())
				return 0;
			return velocityIntegral.evaluateRightHandSideIntegral(cell, shapeFunction1.getVelocityShapeFunction());
		}
	}
}
