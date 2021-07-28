package mixed;

import basic.*;
import tensorproduct.QuadratureRule1D;

import java.util.Objects;

public class MixedCellIntegral<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>,
	PF extends ScalarShapeFunction<CT,FT,ET>, VF extends VectorShapeFunction<CT,FT,ET>> extends CellIntegral<CT
	,MixedShapeFunction<CT,FT,ET,PF,VF>>
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
	
	public static <CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>,
		PF extends ScalarShapeFunction<CT,
		FT,ET>,
		VF extends VectorShapeFunction<CT,FT,ET>> MixedCellIntegral<CT,FT,ET, PF, VF> fromPressureIntegral(CellIntegral<CT
		, PF> pressureIntegral)
	{
		return new MixedCellIntegral<>(pressureIntegral, null);
	}
	
	public static <CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>,
		PF extends ScalarShapeFunction<CT,FT,ET>, VF extends VectorShapeFunction<CT,FT,ET>> MixedCellIntegral<CT,FT,ET, PF,
		VF> fromVelocityIntegral(CellIntegral<CT
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
	                                                  MixedShapeFunction<CT,FT,ET, PF, VF> pressureShapeFunction,
	                                                  MixedShapeFunction<CT,FT,ET, PF, VF> velocityShapeFunction)
	{
		throw new UnsupportedOperationException("needs to be overwritten");
	}
	
	@Override
	public double evaluateCellIntegral(CT cell,
	                                   MixedShapeFunction<CT,FT,ET, PF, VF> shapeFunction1,
	                                   MixedShapeFunction<CT,FT,ET, PF, VF> shapeFunction2)
	{
		if (isPressureIntegral())
		{
			if (shapeFunction1.isVelocity() || shapeFunction2.isVelocity())
				return 0;
			return Objects.requireNonNull(pressureIntegral).evaluateCellIntegral(cell, shapeFunction1.getPressureShapeFunction(),
				shapeFunction2.getPressureShapeFunction());
		} else if (isVelocityIntegral())
		{
			if (shapeFunction1.isPressure() || shapeFunction2.isPressure())
				return 0;
			return Objects.requireNonNull(velocityIntegral).evaluateCellIntegral(cell, shapeFunction1.getVelocityShapeFunction(),
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
