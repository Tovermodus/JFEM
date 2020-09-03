package mixed;

import basic.*;

public class MixedCellIntegral<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
	PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
	PI extends CellIntegral<CT,FT,PF>, VI extends CellIntegral<CT,FT,VF>, PVI extends CellIntegral<CT,FT,ST>>
{
	private final CellIntegral<CT,FT,PF> pressureIntegral;
	private final CellIntegral<CT,FT,VF> velocityIntegral;
	private final CellIntegral<CT,FT,ST> pressureVelocityIntegral;
	
	private MixedCellIntegral(CellIntegral<CT, FT, PF> pressureIntegral,
	                          CellIntegral<CT, FT, VF> velocityIntegral, CellIntegral<CT,FT,ST> pressureVelocityIntegral)
	{
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
		this.pressureVelocityIntegral = pressureVelocityIntegral;
	}
	
	public static<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
		ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
		PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
		PI extends CellIntegral<CT,FT,PF>, VI extends CellIntegral<CT,FT,VF>, PVI extends CellIntegral<CT,FT,
		ST>> MixedCellIntegral<CT,FT,ST,PF,VF,PI,VI,PVI> fromPressureIntegral(CellIntegral<CT,FT
		,PF> pressureIntegral)
	{
		return new MixedCellIntegral<>(pressureIntegral, null, null);
	}
	public static<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
		ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
		PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
		PI extends CellIntegral<CT,FT,PF>, VI extends CellIntegral<CT,FT,VF>, PVI extends CellIntegral<CT,FT,ST>> MixedCellIntegral<CT,FT,ST,PF,
		VF,PI,VI, PVI> fromVelocityIntegral(CellIntegral<CT,FT
		,VF> velocityIntegral)
	{
		return new MixedCellIntegral<>( null, velocityIntegral, null);
	}
	public static<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
		ST extends MixedShapeFunction<CT,FT,ST,PF,VF>,
		PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>,
		PI extends CellIntegral<CT,FT,PF>, VI extends CellIntegral<CT,FT,VF>,PVI extends CellIntegral<CT,FT,ST>> MixedCellIntegral<CT,FT,ST,PF,
		VF,PI,VI,PVI> fromPressureVelocityIntegral(CellIntegral<CT,FT
		,ST> pressureVelocityIntegral)
	{
		return new MixedCellIntegral<>( null, null, pressureVelocityIntegral);
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
	public double evaluateCellIntegral(CT cell,
	                                   ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		if(isPressureIntegral())
			return pressureIntegral.evaluateCellIntegral(cell,shapeFunction1.getPressureFunction(),
				shapeFunction2.getPressureFunction());
		else if(isVelocityIntegral())
			return velocityIntegral.evaluateCellIntegral(cell, shapeFunction1.getVelocityFunction(),
				shapeFunction2.getVelocityFunction());
		else
			return pressureVelocityIntegral.evaluateCellIntegral(cell, shapeFunction1, shapeFunction2);
	}
}
