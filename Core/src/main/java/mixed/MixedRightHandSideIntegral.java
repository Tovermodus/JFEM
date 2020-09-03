package mixed;

import basic.*;

public class MixedRightHandSideIntegral<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
	PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
	PI extends RightHandSideIntegral<CT,FT,PF>, VI extends RightHandSideIntegral<CT,FT,VF>> extends RightHandSideIntegral<CT,FT,ST>
{
	private final RightHandSideIntegral<CT,FT,PF> pressureIntegral;
	private final RightHandSideIntegral<CT,FT,VF> velocityIntegral;
	
	private MixedRightHandSideIntegral(RightHandSideIntegral<CT, FT, PF> pressureIntegral,
	                                   RightHandSideIntegral<CT, FT, VF> velocityIntegral)
	{
		super();
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	public static<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
		ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
		PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
		PI extends RightHandSideIntegral<CT,FT,PF>, VI extends RightHandSideIntegral<CT,FT,VF>> MixedRightHandSideIntegral<CT,FT,ST,PF,VF,PI,VI> fromPressureIntegral(RightHandSideIntegral<CT,FT
		,PF> pressureIntegral)
	{
		return new MixedRightHandSideIntegral<>(pressureIntegral, null);
	}
	public static<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
		ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
		PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
		PI extends RightHandSideIntegral<CT,FT,PF>, VI extends RightHandSideIntegral<CT,FT,VF>> MixedRightHandSideIntegral<CT,FT,ST,PF,
						VF,PI,VI> fromVelocityIntegral(RightHandSideIntegral<CT,FT
		,VF> velocityIntegral)
	{
		return new MixedRightHandSideIntegral<>( null, velocityIntegral);
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
	                                   ST shapeFunction1)
	{
		if(isPressureIntegral())
			return pressureIntegral.evaluateRightHandSideIntegral(cell,shapeFunction1.getPressureFunction());
		else
			return velocityIntegral.evaluateRightHandSideIntegral(cell, shapeFunction1.getVelocityFunction());
	}
}
