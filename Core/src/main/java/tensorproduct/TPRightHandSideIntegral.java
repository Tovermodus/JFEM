package tensorproduct;

import basic.*;
import linalg.CoordinateVector;

public class TPRightHandSideIntegral extends RightHandSideIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>
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
	public double evaluateRightHandSideIntegral(TPCell<TPShapeFunction> cell, TPShapeFunction shapeFunction1)
	{
		if(name.equals(VALUE))
		{
			if(weightIsTensorProduct)
				return TPCellIntegral.integrateTensorProduct(x->shapeFunction1.value(x)*(Double)(rightHandSide.value(x)),cell.cell1Ds);
			else
				return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x)*(Double)(rightHandSide.value(x)),cell.cell1Ds);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
}
