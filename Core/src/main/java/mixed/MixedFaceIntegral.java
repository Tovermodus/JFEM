package mixed;

import basic.Face;
import basic.FaceIntegral;
import basic.ScalarShapeFunction;
import basic.VectorShapeFunction;

public class MixedFaceIntegral<FT extends Face<?, FT>,
	PF extends ScalarShapeFunction<?, FT>, VF extends VectorShapeFunction<?, FT>,
	MF extends ComposeMixedShapeFunction<?, FT, PF, VF>>
	extends FaceIntegral<FT, MF>
{
	private final FaceIntegral<FT, PF> pressureIntegral;
	private final FaceIntegral<FT, VF> velocityIntegral;
	private final FaceIntegral<FT, MF> pressureVelocityIntegral;
	
	private MixedFaceIntegral(final FaceIntegral<FT, PF> pressureIntegral,
	                          final FaceIntegral<FT, VF> velocityIntegral,
	                          final FaceIntegral<FT, MF> pressureVelocityIntegral)
	{
		super(null);
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
		this.pressureVelocityIntegral = pressureVelocityIntegral;
	}
	
	public static <FT extends Face<?, FT>,
		PF extends ScalarShapeFunction<?, FT>, VF extends VectorShapeFunction<?, FT>, MF extends ComposeMixedShapeFunction<?, FT, PF, VF>>
	MixedFaceIntegral<FT, PF, VF, MF>
	fromPressureIntegral(final FaceIntegral<FT, PF> pressureIntegral)
	{
		return new MixedFaceIntegral<>(pressureIntegral, null, null);
	}
	
	public static <FT extends Face<?, FT>,
		PF extends ScalarShapeFunction<?, FT>, VF extends VectorShapeFunction<?, FT>, MF extends ComposeMixedShapeFunction<?, FT, PF, VF>
		> MixedFaceIntegral<FT, PF,
		VF, MF> fromVelocityIntegral(final FaceIntegral<FT
		, VF> velocityIntegral)
	{
		return new MixedFaceIntegral<>(null, velocityIntegral, null);
	}
	
	public static <FT extends Face<?, FT>,
		ST extends ComposeMixedShapeFunction<?, FT, PF, VF>,
		PF extends ScalarShapeFunction<?, FT>, VF extends VectorShapeFunction<?, FT>,
		MF extends ComposeMixedShapeFunction<?, FT, PF, VF>> MixedFaceIntegral<FT, PF, VF, MF> fromPressureVelocityIntegral(
		final FaceIntegral<FT
			, MF> pressureVelocityIntegral)
	{
		return new MixedFaceIntegral<>(null, null, pressureVelocityIntegral);
	}
	
	public boolean isPressureIntegral()
	{
		return this.pressureIntegral != null;
	}
	
	public boolean isVelocityIntegral()
	{
		return this.velocityIntegral != null;
	}
	
	public boolean isPressureVelocityIntegral()
	{
		return this.pressureVelocityIntegral != null;
	}
	
	@Override
	public double evaluateFaceIntegral(final FT face,
	                                   final MF shapeFunction1,
	                                   final MF shapeFunction2)
	{
		if (isPressureIntegral())
			return pressureIntegral.evaluateFaceIntegral(face, shapeFunction1.getPressureFunction(),
			                                             shapeFunction2.getPressureFunction());
		else if (isVelocityIntegral())
			return velocityIntegral.evaluateFaceIntegral(face, shapeFunction1.getVelocityFunction(),
			                                             shapeFunction2.getVelocityFunction());
		else
			return pressureVelocityIntegral.evaluateFaceIntegral(face, shapeFunction1, shapeFunction2);
	}
}
