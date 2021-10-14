package tensorproduct;

import basic.Function;
import basic.RightHandSideIntegral;
import basic.VectorShapeFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TPVectorRightHandSideIntegral<ST extends VectorShapeFunction<TPCell, TPFace>> extends RightHandSideIntegral<TPCell,
	ST>
{
	public static final String VALUE = "Value";
	public static final String H1 = "h1";
	
	public TPVectorRightHandSideIntegral(final Function<?, ?, ?> rightHandSide, final String name)
	{
		super(rightHandSide, name);
		if (name.equals(VALUE) && !(rightHandSide.value(
			new CoordinateVector(rightHandSide.getDomainDimension())) instanceof CoordinateVector))
			throw new IllegalArgumentException();
	}
	
	@Override
	public double evaluateRightHandSideIntegral(final TPCell cell, final ST shapeFunction1)
	{
		if (name.equals(VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(
				x -> shapeFunction1.value(x).inner((CoordinateVector) (rightHandSide.value(x))),
				cell,
				quadratureRule1D);
		}
		if (name.equals(H1))
		{
			//System.out.println(rightHandSide);
//			if (rightHandSide instanceof FunctionOnCells)
//			{
//				final FunctionOnCells<TPCell, ?, ?, ?, ?> rhs
//					= (FunctionOnCells<TPCell, ?, ?, ?, ?>) rightHandSide;
//				return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1
//					                                                .gradientInCell(x, cell)
//					                                                .frobeniusInner(
//						                                                (CoordinateMatrix) rhs.gradientInCell(
//							                                                x, cell)) + shapeFunction1
//					                                                .valueInCell(x, cell)
//					                                                .inner((CoordinateVector) (rhs.valueInCell(
//						                                                x, cell))),
//				                                                cell,
//				                                                quadratureRule1D);
//			}
			return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1
				                                                .gradient(x)
				                                                .frobeniusInner((CoordinateMatrix) rightHandSide.gradient(x)) + shapeFunction1
				                                                .value(x)
				                                                .inner((CoordinateVector) (rightHandSide.value(x))),
			                                                cell,
			                                                quadratureRule1D);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
}
