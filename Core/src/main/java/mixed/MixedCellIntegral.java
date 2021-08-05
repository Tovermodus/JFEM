package mixed;

import basic.*;
import tensorproduct.QuadratureRule1D;

import java.util.Objects;

public class MixedCellIntegral<CT extends Cell<CT,?>,
	PF extends ScalarShapeFunction<CT,?>, VF extends VectorShapeFunction<CT,?>, MF extends MixedShapeFunction<CT,
	?, PF,VF>> extends CellIntegral<CT
	,MF>
{
	private final CellIntegral<CT, PF> pressureIntegral;
	private final CellIntegral<CT, VF> velocityIntegral;
	
	protected MixedCellIntegral(CellIntegral<CT, PF> pressureIntegral,
	                            CellIntegral<CT, VF> velocityIntegral)
	{
		super(null);
		this.pressureIntegral = pressureIntegral;
		this.velocityIntegral = velocityIntegral;
	}
	
	protected MixedCellIntegral(Function<?, ?, ?> weight, String name, QuadratureRule1D quadratureRule1D)
	{
		super(weight, name, quadratureRule1D);
		this.pressureIntegral = null;
		this.velocityIntegral = null;
	}
	
	protected MixedCellIntegral(String name, QuadratureRule1D quadratureRule1D)
	{
		super(name, quadratureRule1D);
		this.pressureIntegral = null;
		this.velocityIntegral = null;
	}
	protected MixedCellIntegral(Function<?, ?, ?> weight, String name)
	{
		this(weight,name,QuadratureRule1D.Gauss5);
	}
	
	protected MixedCellIntegral(String name)
	{
		this(name,QuadratureRule1D.Gauss5);
	}
	
	public static <CT extends Cell<CT, ?>,
		PF extends ScalarShapeFunction<CT,?>,
		VF extends VectorShapeFunction<CT, ?>, MF extends MixedShapeFunction<CT,
		?, PF,VF>> MixedCellIntegral<CT, PF, VF,MF> fromPressureIntegral(CellIntegral<CT
		, PF> pressureIntegral)
	{
		return new MixedCellIntegral<>(pressureIntegral, null);
	}
	
	public static <CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
		PF extends ScalarShapeFunction<CT,FT>, VF extends VectorShapeFunction<CT,FT>,
		MF extends MixedShapeFunction<CT,
	?, PF,VF>> MixedCellIntegral<CT,
		PF,
		VF,MF> fromVelocityIntegral(CellIntegral<CT
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
	                                                  MF pressureShapeFunction,
	                                                  MF velocityShapeFunction)
	{
		throw new UnsupportedOperationException("needs to be overwritten");
	}
	
	@Override
	public double evaluateCellIntegral(CT cell,
	                                   MF shapeFunction1,
	                                   MF shapeFunction2)
	{
		if (isPressureIntegral())
		{
			if (shapeFunction1.hasVelocityFunction() || shapeFunction2.hasVelocityFunction())
				return 0;
			return Objects.requireNonNull(pressureIntegral).evaluateCellIntegral(cell, shapeFunction1.getPressureShapeFunction(),
				shapeFunction2.getPressureShapeFunction());
		} else if (isVelocityIntegral())
		{
			if (shapeFunction1.hasPressureFunction() || shapeFunction2.hasPressureFunction())
				return 0;
			return Objects.requireNonNull(velocityIntegral).evaluateCellIntegral(cell, shapeFunction1.getVelocityShapeFunction(),
				shapeFunction2.getVelocityShapeFunction());
		} else
		{
			if (shapeFunction1.hasPressureFunction() && shapeFunction2.hasVelocityFunction())
				return evaluatePressureVelocityIntegral(cell, shapeFunction2, shapeFunction1);
			else if (shapeFunction2.hasPressureFunction() && shapeFunction1.hasVelocityFunction())
				return evaluatePressureVelocityIntegral(cell, shapeFunction1, shapeFunction2);
			else
				return 0;
		}
	}
}
