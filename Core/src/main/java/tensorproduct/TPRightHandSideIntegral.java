package tensorproduct;

import basic.*;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TPRightHandSideIntegral<ST extends ScalarShapeFunction<TPCell, TPFace>> extends RightHandSideIntegral<TPCell,ST>
{
	public static final String VALUE = "Value";
	public static final String H1 = "H1";
	
	public TPRightHandSideIntegral(Function<?, ?, ?> rightHandSide, String name)
	{
		super(rightHandSide, name);
		if (name.equals(VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(H1) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
	}
	
	@Override
	public double evaluateRightHandSideIntegral(TPCell cell, ST shapeFunction1)
	{
		if (name.equals(VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1.value(x)
					* (Double) (rightHandSide.value(x)),
				cell,
				quadratureRule1D);
		}
		if (name.equals(H1))
		{
			return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1.value(x)
					* (Double) (rightHandSide.value(x))
					+ shapeFunction1.gradient(x).inner((CoordinateVector) (rightHandSide.gradient(x))),
				cell,
				quadratureRule1D);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
}
