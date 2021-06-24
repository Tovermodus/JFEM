package tensorproduct;

import basic.CellIntegral;
import basic.Function;
import basic.ScalarShapeFunction;
import linalg.CoordinateVector;
import linalg.Vector;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class TPCellIntegralViaReferenceCell<ST extends ScalarShapeFunction<TPCell,TPFace,TPEdge,ST>> extends TPCellIntegral<ST>
{
	private double valueOnReferenceCell = Double.NaN;
	public TPCellIntegralViaReferenceCell(Function<?, ?, ?> weight, String name, boolean weightIsTensorProduct)
	{
		super(weight, name, weightIsTensorProduct);
	}
	
	public TPCellIntegralViaReferenceCell(String name)
	{
		super(name);
	}
	@Override
	public double evaluateCellIntegral(TPCell cell, ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		if(name.equals(GRAD_GRAD))
		{
			if(Double.isNaN(valueOnReferenceCell))
				valueOnReferenceCell =
					TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.gr(x).inner(shapeFunction2.gradient(x))*(Double)weight.value(x),
					cell.getReferenceCell());
			return valueOnReferenceCell*;
		}
		return super.evaluateCellIntegral(cell,shapeFunction1,shapeFunction2);
	}
}