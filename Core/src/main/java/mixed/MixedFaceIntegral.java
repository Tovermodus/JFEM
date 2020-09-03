package mixed;

import basic.*;

public class MixedFaceIntegral<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
	PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
	PI extends FaceIntegral<CT,FT,PF>, VI extends FaceIntegral<CT,FT,VF>> extends FaceIntegral<CT,FT,ST>
{
	private final FaceIntegral<CT,FT,PF> pressureIntegral;
	private final FaceIntegral<CT,FT,VF> velocityIntegral;
	private final FaceIntegral<CT,FT,ST> pressureVelocityIntegral;
	
	private MixedFaceIntegral(FaceIntegral<CT, FT, PF> pressureIntegral,
	                          FaceIntegral<CT, FT, VF> velocityIntegral,
	                          FaceIntegral<CT,FT,ST> pressureVelocityIntegral)
	{
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
		this.pressureVelocityIntegral = pressureVelocityIntegral;
	}
	
	public static<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
		ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
		PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
		PI extends FaceIntegral<CT,FT,PF>, VI extends FaceIntegral<CT,FT,VF>> MixedFaceIntegral<CT,FT,ST,PF,VF,PI,VI> fromPressureIntegral(FaceIntegral<CT,FT
		,PF> pressureIntegral)
	{
		return new MixedFaceIntegral<>(pressureIntegral, null, null);
	}
	public static<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
		ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
		PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
		PI extends FaceIntegral<CT,FT,PF>, VI extends FaceIntegral<CT,FT,VF>> MixedFaceIntegral<CT,FT,ST,PF,
		VF,PI,VI> fromVelocityIntegral(FaceIntegral<CT,FT
		,VF> velocityIntegral)
	{
		return new MixedFaceIntegral<>( null, velocityIntegral, null);
	}
	public static<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
		ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
		PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
		PI extends FaceIntegral<CT,FT,PF>, VI extends FaceIntegral<CT,FT,VF>> MixedFaceIntegral<CT,FT,ST,PF,
		VF,PI,VI> fromPressureVelocityIntegral(FaceIntegral<CT,FT
		,ST> pressureVelocityIntegral)
	{
		return new MixedFaceIntegral<>( null, null, pressureVelocityIntegral);
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
	                                   ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		if(isPressureIntegral())
			return pressureIntegral.evaluateFaceIntegral(face,shapeFunction1.getPressureFunction(),
				shapeFunction2.getPressureFunction());
		else if(isVelocityIntegral())
			return velocityIntegral.evaluateFaceIntegral(face, shapeFunction1.getVelocityFunction(),
				shapeFunction2.getVelocityFunction());
		else
			return pressureVelocityIntegral.evaluateFaceIntegral(face, shapeFunction1, shapeFunction2);
	}
}
