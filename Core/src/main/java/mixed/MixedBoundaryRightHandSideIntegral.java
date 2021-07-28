package mixed;

import basic.*;

public class MixedBoundaryRightHandSideIntegral<CT extends Cell<CT,FT,ET>,FT extends Face<CT,FT,ET>,
	ET extends Edge<CT,FT,ET>,
	PF extends ScalarShapeFunction<CT,FT,ET>, VF extends VectorShapeFunction<CT,FT,ET>> extends BoundaryRightHandSideIntegral<FT,MixedShapeFunction<CT,FT,ET,PF,VF>>
{
	private final BoundaryRightHandSideIntegral<FT, PF> pressureIntegral;
	private final BoundaryRightHandSideIntegral<FT, VF> velocityIntegral;
	
	private MixedBoundaryRightHandSideIntegral(BoundaryRightHandSideIntegral<FT, PF> pressureIntegral,
	                                           BoundaryRightHandSideIntegral<FT, VF> velocityIntegral)
	{
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	public static <CT extends Cell<CT, FT, ET>, FT extends Face<CT, FT, ET>, ET extends Edge<CT, FT, ET>,
		PF extends ScalarShapeFunction<CT, FT, ET>, VF extends VectorShapeFunction<CT, FT, ET>
		> MixedBoundaryRightHandSideIntegral<CT, FT, ET, PF,
		VF> fromPressureIntegral(BoundaryRightHandSideIntegral<FT
		, PF> pressureIntegral)
	{
		return new MixedBoundaryRightHandSideIntegral<>(pressureIntegral, null);
	}
	
	public static <CT extends Cell<CT, FT, ET>, FT extends Face<CT, FT, ET>, ET extends Edge<CT, FT, ET>,
		PF extends ScalarShapeFunction<CT, FT, ET>, VF extends VectorShapeFunction<CT, FT, ET>
		> MixedBoundaryRightHandSideIntegral<CT, FT, ET, PF,
		VF> fromVelocityIntegral(BoundaryRightHandSideIntegral<FT
		, VF> velocityIntegral)
	{
		return new MixedBoundaryRightHandSideIntegral<>(null, velocityIntegral);
	}
	
	public boolean isPressureIntegral()
	{
		return this.pressureIntegral != null;
	}
	
	public boolean isVelocityIntegral()
	{
		return this.velocityIntegral != null;
	}
	
	public double evaluateBoundaryRightHandSideIntegral(FT face,
	                                                    MixedShapeFunction<CT, FT, ET, PF, VF> shapeFunction1)
	{
		if (isPressureIntegral())
		{
			if (!shapeFunction1.isPressure())
				return 0;
			return pressureIntegral.evaluateBoundaryRightHandSideIntegral(face,
				shapeFunction1.getPressureShapeFunction());
		} else
		{
			if (!shapeFunction1.isVelocity())
				return 0;
			return velocityIntegral.evaluateBoundaryRightHandSideIntegral(face,
				shapeFunction1.getVelocityShapeFunction());
		}
	}
}
