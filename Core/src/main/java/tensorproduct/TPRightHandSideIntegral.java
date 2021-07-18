package tensorproduct;

import basic.*;
import linalg.CoordinateVector;

public class TPRightHandSideIntegral<ST extends ScalarShapeFunction<TPCell,TPFace,TPEdge>> extends RightHandSideIntegral<TPCell,ST>
{
	public static final String VALUE="Value";
	public static final String H1 = "H1";
	private final boolean weightIsTensorProduct;
	public TPRightHandSideIntegral(Function<?,?,?> rightHandSide, String name, boolean weightIsTensorProduct)
	{
		super(rightHandSide, name);
		if(name.equals(VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(H1) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
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
		if(name.equals(H1))
		{
			return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x)*(Double)(rightHandSide.value(x))
				+ shapeFunction1.gradient(x).inner((CoordinateVector)(rightHandSide.gradient(x))),
				cell.cell1Ds);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
}
