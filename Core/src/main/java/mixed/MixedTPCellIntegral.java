package mixed;

import basic.CellIntegral;
import basic.Function;
import basic.ScalarShapeFunction;
import basic.VectorShapeFunction;
import tensorproduct.ScalarRTFESpace;
import tensorproduct.TPCell;
import tensorproduct.TPCellIntegral;
import tensorproduct.TPFace;

public class MixedTPCellIntegral<PF extends ScalarShapeFunction<TPCell, TPFace, PF>,
	VF extends VectorShapeFunction<TPCell, TPFace, VF>> extends MixedCellIntegral<TPCell,
	TPFace,
	PF, VF>
{
	
	public static final String DIV_VALUE = "DivValue";
	private MixedTPCellIntegral(CellIntegral<TPCell, TPFace, PF> pressureIntegral,
	                            CellIntegral<TPCell, TPFace, VF> velocityIntegral)
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
	protected double evaluatePressureVelocityIntegral(TPCell cell, MixedShapeFunction<TPCell, TPFace, PF, VF> shapeFunction1, MixedShapeFunction<TPCell, TPFace, PF, VF> shapeFunction2)
	{
		if(!isPressureVelocityIntegral())
			throw new IllegalStateException("not a pressure velocity integral");
		if(!shapeFunction1.isVelocity() || !shapeFunction2.isPressure())
			throw new IllegalArgumentException("shapefunctions are not the right type");
		
		if(name.equals(DIV_VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.getVelocityShapeFunction().divergence(x)*shapeFunction2.getPressureFunction().value(x),
				cell.getCell1Ds());
		}
		throw new IllegalStateException("Integral type not supported");
	}
}
