package mixed;

import basic.*;
import com.google.common.base.Stopwatch;

public class MixedCellIntegral<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	PF extends ScalarShapeFunction<CT,FT,PF>, VF extends VectorShapeFunction<CT,FT,VF>> extends CellIntegral<CT,FT,MixedShapeFunction<CT,FT,PF,VF>>
{
	private final CellIntegral<CT, FT, PF> pressureIntegral;
	private final CellIntegral<CT, FT, VF> velocityIntegral;
	
	protected MixedCellIntegral(CellIntegral<CT, FT, PF> pressureIntegral,
	                            CellIntegral<CT, FT, VF> velocityIntegral)
	{
		super();
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	protected MixedCellIntegral(Function<?, ?, ?> weight, String name)
	{
		super(weight, name);
		this.pressureIntegral = null;
		this.velocityIntegral = null;
	}
	
	protected MixedCellIntegral(String name)
	{
		super(name);
		this.pressureIntegral = null;
		this.velocityIntegral = null;
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>> MixedCellIntegral<CT, FT, PF, VF> fromPressureIntegral(CellIntegral<CT, FT
		, PF> pressureIntegral)
	{
		return new MixedCellIntegral<>(pressureIntegral, null);
	}
	
	public static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
		ST extends MixedShapeFunction<CT, FT, PF, VF>,
		PF extends ScalarShapeFunction<CT, FT, PF>, VF extends VectorShapeFunction<CT, FT, VF>> MixedCellIntegral<CT, FT, PF,
		VF> fromVelocityIntegral(CellIntegral<CT, FT
		, VF> velocityIntegral)
	{
		return new MixedCellIntegral<>(null, velocityIntegral);
	}
	
	public boolean isPressureIntegral()
	{
		return pressureIntegral != null;
	}
	
	public boolean isVelocityIntegral()
	{
		return velocityIntegral != null;
	}
	
	public boolean isPressureVelocityIntegral()
	{
		return this.pressureIntegral == null && this.velocityIntegral == null;
	}
	
	protected double evaluatePressureVelocityIntegral(CT cell,
	                                                  MixedShapeFunction<CT, FT, PF, VF> pressureShapeFunction,
	                                                  MixedShapeFunction<CT, FT, PF, VF> velocityShapeFunction)
	{
		throw new UnsupportedOperationException("needs to be overwritten");
	}
	
	@Override
	public double evaluateCellIntegral(CT cell,
	                                   MixedShapeFunction<CT, FT, PF, VF> shapeFunction1,
	                                   MixedShapeFunction<CT, FT, PF, VF> shapeFunction2)
	{
		if (isPressureIntegral())
		{
			if (shapeFunction1.isVelocity() || shapeFunction2.isVelocity())
				return 0;
			return pressureIntegral.evaluateCellIntegral(cell, shapeFunction1.getPressureShapeFunction(),
				shapeFunction2.getPressureShapeFunction());
		} else if (isVelocityIntegral())
		{
			if (shapeFunction1.isPressure() || shapeFunction2.isPressure())
				return 0;
			return velocityIntegral.evaluateCellIntegral(cell, shapeFunction1.getVelocityShapeFunction(),
				shapeFunction2.getVelocityShapeFunction());
		} else
		{
			if (shapeFunction1.isPressure() && shapeFunction2.isVelocity())
				return evaluatePressureVelocityIntegral(cell, shapeFunction2, shapeFunction1);
			else if (shapeFunction2.isPressure() && shapeFunction1.isVelocity())
				return evaluatePressureVelocityIntegral(cell, shapeFunction1, shapeFunction2);
			else
				return 0;
		}
	}
}
