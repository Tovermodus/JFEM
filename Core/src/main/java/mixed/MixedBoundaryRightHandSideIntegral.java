package mixed;

import basic.*;

public class MixedBoundaryRightHandSideIntegral<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>> extends BoundaryRightHandSideIntegral<CT,FT,MixedShapeFunction<CT,FT,PF,VF>>
{
	private final BoundaryRightHandSideIntegral<CT, FT, PF> pressureIntegral;
	private final BoundaryRightHandSideIntegral<CT, FT, VF> velocityIntegral;
	
	private MixedBoundaryRightHandSideIntegral(BoundaryRightHandSideIntegral<CT, FT, PF> pressureIntegral,
	                                           BoundaryRightHandSideIntegral<CT, FT, VF> velocityIntegral)
	{
		super();
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>
		> MixedBoundaryRightHandSideIntegral<CT, FT, PF, VF> fromPressureIntegral(BoundaryRightHandSideIntegral<CT, FT
		, PF> pressureIntegral)
	{
		return new MixedBoundaryRightHandSideIntegral<>(pressureIntegral, null);
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>> MixedBoundaryRightHandSideIntegral<CT, FT, PF,
		VF> fromVelocityIntegral(BoundaryRightHandSideIntegral<CT, FT
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
	                                                    MixedShapeFunction<CT, FT, PF, VF> shapeFunction1)
	{
		if (isPressureIntegral())
			return pressureIntegral.evaluateBoundaryRightHandSideIntegral(face,
				shapeFunction1.getPressureShapeFunction());
		else
			return velocityIntegral.evaluateBoundaryRightHandSideIntegral(face,
				shapeFunction1.getVelocityShapeFunction());
	}
}
