package tensorproduct;

import basic.CellIntegral;
import basic.FunctionOnCells;
import basic.ScalarFunctionOnCells;
import basic.VectorShapeFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TPVectorCellIntegralOnCell<ST extends VectorShapeFunction<TPCell, TPFace>>
	extends CellIntegral<TPCell,
	ST>
{
	public static final String GRAD_GRAD = "GradGrad";
	public static final String VALUE_VALUE = "ValueValue";
	public static final String GRAD_VALUE = "GradValue";
	public static final String VALUE_GRAD = "ValueGrad";
	public static final String SYM_GRAD = "SymGrad";
	public static final String H1 = "H1";
	FunctionOnCells<TPCell, TPFace, ?, ?, ?> weightOnCell;
	
	public TPVectorCellIntegralOnCell(final FunctionOnCells<TPCell, TPFace, ?, ?, ?> weight, final String name)
	{
		this(weight, name, QuadratureRule1D.Gauss5);
	}
	
	public TPVectorCellIntegralOnCell(final double value, final String name)
	{
		this(ScalarFunctionOnCells.constantFunction(value), name);
	}
	
	public TPVectorCellIntegralOnCell(final String name)
	{
		this(name, QuadratureRule1D.Gauss5);
	}
	
	public TPVectorCellIntegralOnCell(final FunctionOnCells<TPCell, TPFace, ?, ?, ?> weight,
	                                  final String name,
	                                  final QuadratureRule1D quadratureRule1D)
	{
		super(weight, name, quadratureRule1D);
		weightOnCell = weight;
		if (name.equals(GRAD_GRAD) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(VALUE_VALUE) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(GRAD_VALUE) && !(weight.defaultValue() instanceof CoordinateVector))
			throw new IllegalArgumentException();
	}
	
	public TPVectorCellIntegralOnCell(final String name, final QuadratureRule1D quadratureRule1D)
	{
		this(ScalarFunctionOnCells.constantFunction(1), name, quadratureRule1D);
	}
	
	@Override
	public double evaluateCellIntegral(final TPCell cell, final ST shapeFunction1, final ST shapeFunction2)
	{
		if (name.equals(SYM_GRAD))
		{
			return TPCellIntegral
				.integrateNonTensorProduct(x ->
				                           {
					                           final CoordinateMatrix grad1 =
						                           shapeFunction1.gradientInCell(x, cell);
					                           final CoordinateMatrix grad2 =
						                           shapeFunction2.gradientInCell(x, cell);
					
					                           return grad1
						                           .add(grad1.transpose())
						                           .frobeniusInner(grad2.add(grad2.transpose()))
						                           * (Double) weightOnCell.valueInCell(x,
						                                                               cell) / 4;
				                           },
				                           cell,
				                           quadratureRule1D);
		}
		if (name.equals(GRAD_GRAD))
		{
			return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1
				                                                .gradientInCell(x, cell)
				                                                .frobeniusInner(shapeFunction2.gradientInCell(x, cell))
				                                                * (Double) weightOnCell.valueInCell(x, cell),
			                                                cell,
			                                                quadratureRule1D);
		}
		if (name.equals(H1))
		{
			return TPCellIntegral.integrateNonTensorProduct(x -> (shapeFunction1
				                                                      .gradientInCell(x, cell)
				                                                      .frobeniusInner(
					                                                      shapeFunction2.gradientInCell(
						                                                      x, cell))
				                                                      + shapeFunction1
				                                                .valueInCell(x, cell)
				                                                .inner(shapeFunction2.valueInCell(x, cell)))
				                                                * (Double) weightOnCell.valueInCell(x, cell),
			                                                cell,
			                                                quadratureRule1D);
		}
		if (name.equals(VALUE_VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1
				                                                .valueInCell(x, cell)
				                                                .inner(shapeFunction2.valueInCell(x, cell))
				                                                * (Double) weightOnCell.valueInCell(x, cell),
			                                                cell,
			                                                quadratureRule1D);
		}
		if (name.equals(GRAD_VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(x ->
				                                                shapeFunction1
					                                                .gradientInCell(x, cell)
					                                                .transpose()
					                                                .mvMul(shapeFunction2.valueInCell(
						                                                x,
						                                                cell))
					                                                .inner((CoordinateVector) weightOnCell
						                                                .valueInCell(x, cell)),
			                                                cell,
			                                                quadratureRule1D);
		}
		if (name.equals(VALUE_GRAD))
		{
			return TPCellIntegral.integrateNonTensorProduct(x ->
				                                                shapeFunction2
					                                                .gradientInCell(x, cell)
					                                                .mvMul(shapeFunction1.valueInCell(
						                                                x,
						                                                cell))
					                                                .inner((CoordinateVector) weightOnCell
						                                                .valueInCell(x, cell)),
			                                                cell,
			                                                quadratureRule1D);
		}
		throw new UnsupportedOperationException("unknown integral name");
	}
}
