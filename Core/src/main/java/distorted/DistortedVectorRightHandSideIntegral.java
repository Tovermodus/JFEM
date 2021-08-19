package distorted;

import basic.RightHandSideIntegral;
import basic.VectorFunction;
import distorted.geometry.DistortedCell;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public class DistortedVectorRightHandSideIntegral extends RightHandSideIntegral<DistortedCell, DistortedVectorShapeFunction>
{
	public static final String H1 = "H1";
	
	public DistortedVectorRightHandSideIntegral(final VectorFunction rightHandSide, final String name)
	{
		super(rightHandSide, name);
	}
	
	@Override
	public double evaluateRightHandSideIntegral(final DistortedCell cell, final DistortedVectorShapeFunction shapeFunction1)
	{
		if (name.equals(H1))
		{
			return DistortedCellIntegral
				.integrateOnReferenceCell
					(x ->
					 {
						 final CoordinateVector xOnGrid = cell.transformFromReferenceCell(x);
						 return shapeFunction1.valueOnReferenceCell(x, cell).inner(
							 (CoordinateVector) (rightHandSide.value(
								 cell.transformFromReferenceCell(x))))
							 + shapeFunction1.gradientOnReferenceCell(x, cell)
							                 .frobeniusInner(
								                 (CoordinateMatrix) (rightHandSide.gradient(
									                 xOnGrid)));
					 },
					 cell,
					 quadratureRule1D);
		} else throw new IllegalArgumentException("Name unknown");
	}
}
