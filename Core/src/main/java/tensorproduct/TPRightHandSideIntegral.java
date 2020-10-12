package tensorproduct;

import basic.*;
import linalg.CoordinateVector;

public class TPRightHandSideIntegral<ST extends ScalarShapeFunction<TPCell,TPFace,TPEdge,ST>> extends RightHandSideIntegral<TPCell,ST>
{
	public static final String VALUE="Value";
	private final boolean weightIsTensorProduct;
	public TPRightHandSideIntegral(Function<?,?,?> rightHandSide, String name, boolean weightIsTensorProduct)
	{
		super(rightHandSide, name);
		if(name.equals(VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		this.weightIsTensorProduct = weightIsTensorProduct;
	}
	
	@Override
	public double evaluateRightHandSideIntegral(TPCell cell, ST shapeFunction1)
	{
		if(name.equals(VALUE))
		{
				return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x)*(Double)(rightHandSide.value(x)),cell.cell1Ds);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
}
