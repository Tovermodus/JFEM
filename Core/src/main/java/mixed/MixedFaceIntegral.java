package mixed;

import basic.*;

public class MixedFaceIntegral<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>,
	PF extends ScalarShapeFunction<CT,FT,ET>, VF extends VectorShapeFunction<CT,FT,ET>> extends FaceIntegral<FT,MixedShapeFunction<CT,FT,ET,PF,VF>>
{
	private final FaceIntegral<FT,PF> pressureIntegral;
	private final FaceIntegral<FT,VF> velocityIntegral;
	private final FaceIntegral<FT, MixedShapeFunction<CT, FT,ET, PF, VF>> pressureVelocityIntegral;
	
	private MixedFaceIntegral(FaceIntegral<FT, PF> pressureIntegral,
	                          FaceIntegral<FT, VF> velocityIntegral,
	                          FaceIntegral<FT, MixedShapeFunction<CT, FT, ET,PF, VF>> pressureVelocityIntegral)
	{
		super(null);
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
		this.pressureVelocityIntegral = pressureVelocityIntegral;
	}
	
	public static <CT extends Cell<CT, FT,ET>, FT extends Face<CT, FT,ET>, ET extends Edge<CT,FT,ET>,
		PF extends ScalarShapeFunction<CT, FT, ET>, VF extends VectorShapeFunction<CT, FT,ET>> MixedFaceIntegral<CT, FT,ET, PF, VF> fromPressureIntegral(FaceIntegral<FT,PF> pressureIntegral)
	{
		return new MixedFaceIntegral<>(pressureIntegral, null, null);
	}
	
	public static <CT extends Cell<CT, FT,ET>, FT extends Face<CT, FT,ET>, ET extends Edge<CT,FT,ET>,
		ST extends MixedShapeFunction<CT, FT, ET,PF, VF>,
		PF extends ScalarShapeFunction<CT, FT, ET>, VF extends VectorShapeFunction<CT, FT,ET>> MixedFaceIntegral<CT,FT,ET, PF,
		VF> fromVelocityIntegral(FaceIntegral<FT
		, VF> velocityIntegral)
	{
		return new MixedFaceIntegral<>(null, velocityIntegral, null);
	}
	
	public static <CT extends Cell<CT, FT,ET>, FT extends Face<CT, FT,ET>, ET extends Edge<CT,FT,ET>,
		ST extends MixedShapeFunction<CT, FT, ET,PF, VF>,
		PF extends ScalarShapeFunction<CT, FT, ET>, VF extends VectorShapeFunction<CT, FT,ET>> MixedFaceIntegral<CT,FT,ET, PF,
		VF> fromPressureVelocityIntegral(FaceIntegral<FT
		, MixedShapeFunction<CT, FT,ET, PF, VF>> pressureVelocityIntegral)
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
	                                   MixedShapeFunction<CT, FT,ET, PF, VF> shapeFunction1,
	                                   MixedShapeFunction<CT, FT, ET,PF, VF> shapeFunction2)
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
