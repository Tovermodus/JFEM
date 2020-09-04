package mixed;

import basic.*;

public class MixedFaceIntegral<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>> extends FaceIntegral<CT,FT,
	MixedShapeFunction<CT,FT,PF,VF>>
{
	private final FaceIntegral<CT, FT, PF> pressureIntegral;
	private final FaceIntegral<CT, FT, VF> velocityIntegral;
	private final FaceIntegral<CT, FT, MixedShapeFunction<CT, FT, PF, VF>> pressureVelocityIntegral;
	
	private MixedFaceIntegral(FaceIntegral<CT, FT, PF> pressureIntegral,
	                          FaceIntegral<CT, FT, VF> velocityIntegral,
	                          FaceIntegral<CT, FT, MixedShapeFunction<CT, FT, PF, VF>> pressureVelocityIntegral)
	{
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
		this.pressureVelocityIntegral = pressureVelocityIntegral;
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>> MixedFaceIntegral<CT, FT, PF, VF> fromPressureIntegral(FaceIntegral<CT, FT
		, PF> pressureIntegral)
	{
		return new MixedFaceIntegral<>(pressureIntegral, null, null);
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		ST extends MixedShapeFunction<CT, FT, PF, VF>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>> MixedFaceIntegral<CT, FT, PF,
		VF> fromVelocityIntegral(FaceIntegral<CT, FT
		, VF> velocityIntegral)
	{
		return new MixedFaceIntegral<>(null, velocityIntegral, null);
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		ST extends MixedShapeFunction<CT, FT, PF, VF>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>
		> MixedFaceIntegral<CT, FT, PF,
		VF> fromPressureVelocityIntegral(FaceIntegral<CT, FT
		, MixedShapeFunction<CT, FT, PF, VF>> pressureVelocityIntegral)
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
	public double evaluateFaceIntegral(FT face,
	                                   MixedShapeFunction<CT, FT, PF, VF> shapeFunction1,
	                                   MixedShapeFunction<CT, FT, PF, VF> shapeFunction2)
	{
		if (isPressureIntegral())
			return pressureIntegral.evaluateFaceIntegral(face, shapeFunction1.getPressureShapeFunction(),
				shapeFunction2.getPressureShapeFunction());
		else if (isVelocityIntegral())
			return velocityIntegral.evaluateFaceIntegral(face, shapeFunction1.getVelocityShapeFunction(),
				shapeFunction2.getVelocityShapeFunction());
		else
			return pressureVelocityIntegral.evaluateFaceIntegral(face, shapeFunction1, shapeFunction2);
	}
}
