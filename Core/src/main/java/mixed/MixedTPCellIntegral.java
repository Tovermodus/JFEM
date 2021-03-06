package mixed;

import basic.CellIntegral;
import basic.Function;
import basic.ScalarShapeFunction;
import basic.VectorShapeFunction;
import tensorproduct.*;

public class MixedTPCellIntegral<PF extends ScalarShapeFunction<TPCell, TPFace, TPEdge>,
	VF extends VectorShapeFunction<TPCell, TPFace,TPEdge>> extends MixedCellIntegral<TPCell, TPFace, TPEdge,PF, VF>
{
	
	public static final String DIV_VALUE = "DivValue";
	public static final String VALUE_GRAD = "ValueGrad";
	private MixedTPCellIntegral(CellIntegral<TPCell, PF> pressureIntegral,
	                            CellIntegral<TPCell, VF> velocityIntegral)
	{
		super(pressureIntegral, velocityIntegral);
	}
	
	public MixedTPCellIntegral(Function<?, ?, ?> weight, String name)
	{
		super(weight, name);
	}
	
	public MixedTPCellIntegral(String name)
	{
		super(name);
	}
	
	@Override
	protected double evaluatePressureVelocityIntegral(TPCell cell,
	                                                  MixedShapeFunction<TPCell, TPFace, TPEdge,PF, VF> shapeFunction1
		, MixedShapeFunction<TPCell, TPFace,TPEdge, PF, VF> shapeFunction2)
	{
		if(!isPressureVelocityIntegral())
			throw new IllegalStateException("not a pressure velocity integral");
		if(!shapeFunction1.isVelocity() || !shapeFunction2.isPressure())
			throw new IllegalArgumentException("shapefunctions are not the right type");
		
		if(name.equals(DIV_VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.getVelocityShapeFunction().divergence(x)*shapeFunction2.getPressureFunction().value(x)*(Double)weight.value(x),
				cell.getCell1Ds());
		}
		if(name.equals(VALUE_GRAD))
		{
			return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.getVelocityShapeFunction().value(x).inner(shapeFunction2.getPressureFunction().gradient(x))*(Double)weight.value(x),
				cell.getCell1Ds());
		}
		throw new IllegalStateException("Integral type not supported");
	}
}
