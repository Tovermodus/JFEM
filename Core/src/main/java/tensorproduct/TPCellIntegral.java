package tensorproduct;

import basic.CellIntegral;
import basic.Function;
import basic.ScalarFunction;
import basic.ScalarShapeFunction;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.Vector;
import tensorproduct.geometry.Cell1D;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class TPCellIntegral<ST extends ScalarShapeFunction<TPCell, TPFace>>
	extends CellIntegral<TPCell, ST>
{
	public static final String GRAD_GRAD = "GradGrad";
	public static final String VALUE_VALUE = "ValueValue";
	public static final String GRAD_VALUE = "GradValue";
	public static final String VALUE_GRAD = "ValueGrad";
	public static final String H1 = "H1";
	
	public TPCellIntegral(final double weight, final String name, final QuadratureRule1D quadratureRule1D)
	{
		this(ScalarFunction.constantFunction(weight), name, quadratureRule1D);
	}
	
	public TPCellIntegral(final double weight, final String name)
	{
		this(ScalarFunction.constantFunction(weight), name);
	}
	
	public TPCellIntegral(final Function<?, ?, ?> weight, final String name)
	{
		this(weight, name, QuadratureRule1D.Gauss5);
	}
	
	public TPCellIntegral(final Function<?, ?, ?> weight,
	                      final String name,
	                      final QuadratureRule1D quadratureRule1D)
	{
		super(weight, name, quadratureRule1D);
		if (name.equals(GRAD_GRAD) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(H1) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(VALUE_VALUE) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(GRAD_VALUE) && !(weight.defaultValue() instanceof Vector))
			throw new IllegalArgumentException();
		if (name.equals(VALUE_GRAD) && !(weight.defaultValue() instanceof Vector))
			throw new IllegalArgumentException();
	}
	
	public TPCellIntegral(final String name)
	{
		this(name, QuadratureRule1D.Gauss5);
	}
	
	public TPCellIntegral(final String name, final QuadratureRule1D quadratureRule1D)
	{
		super(name, quadratureRule1D);
		if (name.equals(GRAD_GRAD) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(H1) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(VALUE_VALUE) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(GRAD_VALUE) && !(weight.defaultValue() instanceof Vector))
			throw new IllegalArgumentException();
		if (name.equals(VALUE_GRAD) && !(weight.defaultValue() instanceof Vector))
			throw new IllegalArgumentException();
	}
	
	public static double integrateNonTensorProduct(final ToDoubleFunction<CoordinateVector> eval,
	                                               final List<Cell1D> cells,
	                                               final QuadratureRule1D quadratureRule)
	{
		double ret = 0;
		double val;
		final CoordinateVector quadraturePoint = new CoordinateVector(cells.size());
		final double[][][] pointsWeights = new double[cells.size()][2][quadratureRule.length()];
		for (int j = 0; j < cells.size(); j++)
			pointsWeights[j] = cells.get(j)
			                        .distributeQuadrature(quadratureRule);
		final IntCoordinates quadraturePointSize = IntCoordinates.repeat(quadratureRule.length(), cells.size());
		for (final IntCoordinates c : quadraturePointSize.range())
		{
			for (int i = 0; i < c.getDimension(); i++)
			{
				quadraturePoint.set(pointsWeights[i][0][c.get(i)], i);
			}
			val = eval.applyAsDouble(quadraturePoint);
			for (int i = 0; i < cells.size(); i++)
				val *= pointsWeights[i][1][c.get(i)];
			ret += val;
		}
		return ret;
	}
	
	public static double integrateNonTensorProduct(final ToDoubleFunction<CoordinateVector> eval, final TPCell cell,
	                                               final QuadratureRule1D quadratureRule)
	{
		final List<Cell1D> cells = cell.getComponentCells();
		return integrateNonTensorProduct(eval, cells, quadratureRule);
	}
	
	@Override
	public double evaluateCellIntegral(final TPCell cell, final ST shapeFunction1,
	                                   final ST shapeFunction2)
	{
		if (name.equals(GRAD_GRAD))
		{
			return integrateNonTensorProduct(x -> shapeFunction1
				                                 .gradientInCell(x, cell)
				                                 .inner(shapeFunction2.gradientInCell(x, cell))
				                                 * (Double) weight.value(x),
			                                 cell,
			                                 quadratureRule1D);
		}
		if (name.equals(H1))
		{
			return integrateNonTensorProduct(x -> (shapeFunction1.value(x) * shapeFunction2.value(x)
				                                       + shapeFunction1.gradient(x)
				                                                       .inner(
					                                                       shapeFunction2.gradient(x)))
				                                 * (Double) weight.value(x),
			                                 cell,
			                                 quadratureRule1D);
		}
		if (name.equals(VALUE_VALUE))
		{
			return integrateNonTensorProduct(x -> shapeFunction1.value(x)
				                                 * shapeFunction2.value(x)
				                                 * (Double) weight.value(x),
			                                 cell,
			                                 quadratureRule1D);
		}
		if (name.equals(GRAD_VALUE))
		{
			return integrateNonTensorProduct(x -> shapeFunction1.gradient(x)
			                                                    .inner((Vector) weight.value(x))
				                                 * shapeFunction2.value(x),
			                                 cell,
			                                 quadratureRule1D);
		}
		if (name.equals(VALUE_GRAD))
		{
			return integrateNonTensorProduct(x -> shapeFunction2.gradient(x)
			                                                    .inner((Vector) weight.value(x))
				                                 * shapeFunction1.value(x),
			                                 cell,
			                                 quadratureRule1D);
		}
		throw new UnsupportedOperationException("unknown integral name");
	}
}
