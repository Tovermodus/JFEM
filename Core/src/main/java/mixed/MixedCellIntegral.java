package mixed;

import basic.*;

public class MixedCellIntegral<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>> extends CellIntegral<CT,FT,MixedShapeFunction<CT,FT,PF,VF>>
{
	private final CellIntegral<CT, FT, PF> pressureIntegral;
	private final CellIntegral<CT, FT, VF> velocityIntegral;
	private final CellIntegral<CT, FT, MixedShapeFunction<CT, FT, PF, VF>> pressureVelocityIntegral;
	
	private MixedCellIntegral(CellIntegral<CT, FT, PF> pressureIntegral,
	                          CellIntegral<CT, FT, VF> velocityIntegral, CellIntegral<CT, FT, MixedShapeFunction<CT, FT, PF, VF>> pressureVelocityIntegral)
	{
		super();
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
		this.pressureVelocityIntegral = pressureVelocityIntegral;
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>> MixedCellIntegral<CT, FT, PF, VF> fromPressureIntegral(CellIntegral<CT, FT
		, PF> pressureIntegral)
	{
		return new MixedCellIntegral<>(pressureIntegral, null, null);
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		ST extends MixedShapeFunction<CT, FT, PF, VF>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>> MixedCellIntegral<CT, FT, PF,
		VF> fromVelocityIntegral(CellIntegral<CT, FT
		, VF> velocityIntegral)
	{
		return new MixedCellIntegral<>(null, velocityIntegral, null);
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>> MixedCellIntegral<CT, FT, PF,
		VF> fromPressureVelocityIntegral(CellIntegral<CT, FT
		, MixedShapeFunction<CT, FT, PF, VF>> pressureVelocityIntegral)
	{
		return new MixedCellIntegral<>(null, null, pressureVelocityIntegral);
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
	public double evaluateCellIntegral(CT cell,
	                                   MixedShapeFunction<CT, FT, PF, VF> shapeFunction1,
	                                   MixedShapeFunction<CT, FT, PF, VF> shapeFunction2)
	{
		if (isPressureIntegral())
			return pressureIntegral.evaluateCellIntegral(cell, shapeFunction1.getPressureShapeFunction(),
				shapeFunction2.getPressureShapeFunction());
		else if (isVelocityIntegral())
			return velocityIntegral.evaluateCellIntegral(cell, shapeFunction1.getVelocityShapeFunction(),
				shapeFunction2.getVelocityShapeFunction());
		else
			return pressureVelocityIntegral.evaluateCellIntegral(cell, shapeFunction1, shapeFunction2);
	}
}
