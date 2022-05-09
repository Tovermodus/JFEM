package tensorproduct;

import basic.FunctionOnCells;
import basic.RightHandSideIntegral;
import basic.VectorShapeFunction;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TPVectorRightHandSideIntegralOnCell<ST extends VectorShapeFunction<TPCell, TPFace>>
	extends RightHandSideIntegral<TPCell,
	ST>
{
	public static final String VALUE = "Value";
	public static final String H1 = "h1";
	FunctionOnCells<TPCell, TPFace, ?, ?, ?> rhs;
	
	public TPVectorRightHandSideIntegralOnCell(final FunctionOnCells<TPCell, TPFace, ?, ?, ?> rightHandSide,
	                                           final String name)
	{
		super(rightHandSide, name);
		rhs = rightHandSide;
		if (name.equals(VALUE) && !(rightHandSide.defaultValue() instanceof CoordinateVector))
			throw new IllegalArgumentException();
	}
	
	@Override
	public double evaluateRightHandSideIntegral(final TPCell cell, final ST shapeFunction1)
	{
		if (name.equals(VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(
				x ->
				{
//					System.out.println(shapeFunction1.valueInCell(x,
//					                                              cell) + " " + rhs.valueInCell(x
//						, cell));
					return shapeFunction1.valueInCell(x, cell)
					                     .inner((CoordinateVector) rhs.valueInCell(x,
					                                                               cell));
				},
				cell,
				quadratureRule1D);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
}
