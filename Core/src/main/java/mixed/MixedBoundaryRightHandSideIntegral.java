package mixed;

import basic.*;

public class MixedBoundaryRightHandSideIntegral<FT extends Face<?,FT>,
	PF extends ScalarShapeFunction<?,FT>, VF extends VectorShapeFunction<?,
	FT>, MF extends MixedShapeFunction<?,FT,PF,VF>> extends BoundaryRightHandSideIntegral<FT, MF>
{
	private final BoundaryRightHandSideIntegral<FT, PF> pressureIntegral;
	private final BoundaryRightHandSideIntegral<FT, VF> velocityIntegral;
	
	private MixedBoundaryRightHandSideIntegral(BoundaryRightHandSideIntegral<FT, PF> pressureIntegral,
	                                           BoundaryRightHandSideIntegral<FT, VF> velocityIntegral)
	{
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	public static <FT extends Face<?, FT>,
		PF extends ScalarShapeFunction<?, FT>, VF extends VectorShapeFunction<?, FT>, MF extends MixedShapeFunction<?,FT,PF,VF>
		> MixedBoundaryRightHandSideIntegral<FT, PF,
		VF,MF> fromPressureIntegral(BoundaryRightHandSideIntegral<FT
		, PF> pressureIntegral)
	{
		return new MixedBoundaryRightHandSideIntegral<>(pressureIntegral, null);
	}
	
	public static <FT extends Face<?, FT>,
		PF extends ScalarShapeFunction<?, FT>, VF extends VectorShapeFunction<?, FT>, MF extends MixedShapeFunction<?,FT,PF,VF>
		> MixedBoundaryRightHandSideIntegral<FT, PF,
		VF, MF> fromVelocityIntegral(BoundaryRightHandSideIntegral<FT
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
	                                                    MF shapeFunction1)
	{
		if (isPressureIntegral())
		{
			if (!shapeFunction1.hasPressureFunction())
				return 0;
			return pressureIntegral.evaluateBoundaryRightHandSideIntegral(face,
				shapeFunction1.getPressureShapeFunction());
		} else
		{
			if (!shapeFunction1.hasVelocityFunction())
				return 0;
			return velocityIntegral.evaluateBoundaryRightHandSideIntegral(face,
				shapeFunction1.getVelocityShapeFunction());
		}
	}
}
