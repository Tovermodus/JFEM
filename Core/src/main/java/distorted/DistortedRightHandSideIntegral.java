package distorted;

import basic.RightHandSideIntegral;
import basic.ScalarFunction;
import distorted.geometry.DistortedCell;
import linalg.CoordinateVector;

public class DistortedRightHandSideIntegral extends RightHandSideIntegral<DistortedCell, DistortedShapeFunction>
{
	public static final String VALUE = "VALUE";
	public static final String H1 = "H1";
	
	public DistortedRightHandSideIntegral(final ScalarFunction rightHandSide, final String name)
	{
		super(rightHandSide, name);
	}
	
	@Override
	public double evaluateRightHandSideIntegral(final DistortedCell cell, final DistortedShapeFunction shapeFunction1)
	{
		if (name.equals(VALUE))
		{
			return DistortedCellIntegral
				.integrateOnReferenceCell
					(x -> shapeFunction1
						 .valueOnReferenceCell(x, cell)
						 * (Double) rightHandSide.value(cell.transformFromReferenceCell(x)),
					 cell,
					 quadratureRule1D);
		}
		if (name.equals(H1))
		{
			return DistortedCellIntegral
				.integrateOnReferenceCell
					(x ->
					 {
						 CoordinateVector xOnGrid = cell.transformFromReferenceCell(x);
						 return shapeFunction1.valueOnReferenceCell(x, cell)
							 * (Double) (rightHandSide.value(
							 cell.transformFromReferenceCell(x)))
							 + shapeFunction1.gradientOnReferenceCell(x, cell)
							                 .inner((CoordinateVector) (rightHandSide.gradient(
								                 xOnGrid)));
					 },
					 cell,
					 quadratureRule1D);
		} else throw new IllegalArgumentException("Name unknown");
	}
}
