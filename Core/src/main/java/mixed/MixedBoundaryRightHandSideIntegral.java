package mixed;

import basic.BoundaryRightHandSideIntegral;
import basic.Face;
import basic.ScalarShapeFunction;
import basic.VectorShapeFunction;

public class MixedBoundaryRightHandSideIntegral<FT extends Face<?, FT>,
	PF extends ScalarShapeFunction<?, FT>, VF extends VectorShapeFunction<?,
	FT>, MF extends ComposeMixedShapeFunction<?, FT, PF, VF>>
	extends BoundaryRightHandSideIntegral<FT, MF>
{
	private final BoundaryRightHandSideIntegral<FT, PF> pressureIntegral;
	private final BoundaryRightHandSideIntegral<FT, VF> velocityIntegral;
	
	private MixedBoundaryRightHandSideIntegral(final BoundaryRightHandSideIntegral<FT, PF> pressureIntegral,
	                                           final BoundaryRightHandSideIntegral<FT, VF> velocityIntegral)
	{
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	public static <FT extends Face<?, FT>,
		PF extends ScalarShapeFunction<?, FT>, VF extends VectorShapeFunction<?, FT>, MF extends ComposeMixedShapeFunction<?, FT, PF, VF>
		> MixedBoundaryRightHandSideIntegral<FT, PF,
		VF, MF> fromPressureIntegral(final BoundaryRightHandSideIntegral<FT
		, PF> pressureIntegral)
	{
		return new MixedBoundaryRightHandSideIntegral<>(pressureIntegral, null);
	}
	
	public static <FT extends Face<?, FT>,
		PF extends ScalarShapeFunction<?, FT>, VF extends VectorShapeFunction<?, FT>, MF extends ComposeMixedShapeFunction<?, FT, PF, VF>
		> MixedBoundaryRightHandSideIntegral<FT, PF,
		VF, MF> fromVelocityIntegral(final BoundaryRightHandSideIntegral<FT
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
	
	@Override
	public double evaluateBoundaryRightHandSideIntegral(final FT face,
	                                                    final MF shapeFunction1)
	{
		if (isPressureIntegral())
		{
			if (!shapeFunction1.hasPressureFunction())
				return 0;
			return pressureIntegral.evaluateBoundaryRightHandSideIntegral(face,
			                                                              shapeFunction1.getPressureFunction());
		} else
		{
			if (!shapeFunction1.hasVelocityFunction())
				return 0;
			return velocityIntegral.evaluateBoundaryRightHandSideIntegral(face,
			                                                              shapeFunction1.getVelocityFunction());
		}
	}
}
