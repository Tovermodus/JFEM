package tensorproduct;

import basic.Function;
import basic.RightHandSideIntegral;
import basic.VectorShapeFunction;
import linalg.CoordinateVector;

public class TPVectorRightHandSideIntegral<ST extends VectorShapeFunction<TPCell,TPFace,ST>> extends RightHandSideIntegral<TPCell,
	TPFace,ST>
{
	public static final String VALUE="Value";
	public TPVectorRightHandSideIntegral(Function<?,?,?> rightHandSide, String name, boolean weightIsTensorProduct)
	{
		super(rightHandSide, name);
		if(name.equals(VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof CoordinateVector))
			throw new IllegalArgumentException();
	}
	
	@Override
	public double evaluateRightHandSideIntegral(TPCell cell, ST shapeFunction1)
	{
		if(name.equals(VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x).inner((CoordinateVector)(rightHandSide.value(x))),cell.cell1Ds);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
}
