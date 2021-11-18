package mixed;

import basic.CellIntegral;
import basic.Function;
import basic.ScalarShapeFunction;
import basic.VectorShapeFunction;
import tensorproduct.TPCellIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class MixedTPCellIntegral<PF extends ScalarShapeFunction<TPCell, TPFace>,
	VF extends VectorShapeFunction<TPCell, TPFace>, MF extends ComposeMixedShapeFunction<TPCell, TPFace, PF, VF>>
	extends MixedCellIntegral<TPCell, PF, VF, MF>
{
	
	public static final String DIV_VALUE = "DivValue";
	public static final String VALUE_GRAD = "ValueGrad";
	
	private MixedTPCellIntegral(final CellIntegral<TPCell, PF> pressureIntegral,
	                            final CellIntegral<TPCell, VF> velocityIntegral)
	{
		super(pressureIntegral, velocityIntegral);
	}
	
	public MixedTPCellIntegral(final Function<?, ?, ?> weight, final String name)
	{
		super(weight, name);
	}
	
	public MixedTPCellIntegral(final String name)
	{
		super(name);
	}
	
	@Override
	protected double evaluatePressureVelocityIntegral(final TPCell cell,
	                                                  final MF shapeFunction1
		, final MF shapeFunction2)
	{
		if (!isPressureVelocityIntegral())
			throw new IllegalStateException("not a pressure velocity integral");
		if (!shapeFunction1.hasVelocityFunction() || !shapeFunction2.hasPressureFunction())
			throw new IllegalArgumentException("shapefunctions are not the right type");
		
		if (name.equals(DIV_VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1.getVelocityFunction()
			                                                                   .divergence(x)
				                                                * shapeFunction2.getPressureFunction()
				                                                                .value(x)
				                                                * (Double) weight.value(x),
			                                                cell,
			                                                quadratureRule1D);
		}
		if (name.equals(VALUE_GRAD))
		{
			return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1.getVelocityFunction()
			                                                                   .value(x)
			                                                                   .inner(shapeFunction2.getPressureFunction()
			                                                                                        .gradient(
				                                                                                        x))
				                                                * (Double) weight.value(x),
			                                                cell,
			                                                quadratureRule1D);
		}
		throw new IllegalStateException("Integral type not supported");
	}
}
