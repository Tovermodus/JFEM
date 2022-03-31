package tensorproduct;

import basic.CellIntegral;
import basic.Function;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import tensorproduct.geometry.Cell1D;
import tensorproduct.geometry.TPCell;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class TPCellIntegralOnReferenceCell
	extends CellIntegral<TPCell, CTPShapeFunction>
{
	public static final String GRAD_GRAD = "GradGrad";
	final TPCell referenceCell;
	CoordinateVector[] points;
	double[] weights;
	
	public TPCellIntegralOnReferenceCell(final Function<?, ?, ?> weight, final String name, final int dimension)
	{
		super(weight, name, QuadratureRule1D.Gauss5);
		referenceCell = TPCell.unitHyperCube(dimension);
		final IntCoordinates quadraturePointSize = IntCoordinates.repeat(QuadratureRule1D.Gauss5.length(),
		                                                                 dimension);
		points = new CoordinateVector[quadraturePointSize.size()];
		weights = new double[quadraturePointSize.size()];
		int index = 0;
		for (final IntCoordinates c : quadraturePointSize.range())
		{
			points[index] = new CoordinateVector(dimension);
			weights[index] = 1;
			index++;
		}
		for (int i = 0; i < dimension; i++)
		{
			index = 0;
			final double[][] pointsWeigths =
				referenceCell.getComponentCell(i)
				             .distributeQuadrature(QuadratureRule1D.Gauss5);
			for (final IntCoordinates c : quadraturePointSize.range())
			{
				points[index].set(pointsWeigths[0][c.get(i)], i);
				weights[index] *= pointsWeigths[1][c.get(i)];
				index++;
			}
		}
	}
	
	public double integrateOnReferenceCell(final ToDoubleFunction<CoordinateVector> eval, final TPCell cell)
	{
		double ret = 0;
		double val;
		final int d = referenceCell.getDimension();
		double scaling = 1;
		for (int j = 0; j < d; j++)
		{
			scaling *= cell.size[j];
		}
		for (int i = 0; i < points.length; i++)
		{
//			System.out.println();
//			System.out.println(i);
//			System.out.println("point " + points[i]);
//			System.out.println("value " + eval.applyAsDouble(points[i]));
//			System.out.println("weight " + weights[i]);
//			System.out.println("scaling " + scaling);
			
			val = eval.applyAsDouble(points[i])
				* scaling * weights[i];
			ret += val;
		}
		//System.out.println(ret);
		return ret;
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
			System.out.println();
			System.out.println(c);
			System.out.println("point " + quadraturePoint);
			val = eval.applyAsDouble(quadraturePoint);
			System.out.println("value " + val);
			double w = 1;
			for (int i = 0; i < cells.size(); i++)
				w *= pointsWeights[i][1][c.get(i)];
			System.out.println("weight " + w);
			ret += val * w;
		}
		System.out.println(ret);
		return ret;
	}
	
	@Override
	public double evaluateCellIntegral(final TPCell cell, final CTPShapeFunction shapeFunction1,
	                                   final CTPShapeFunction shapeFunction2)
	{
		if (name.equals(GRAD_GRAD))
		{
//			integrateNonTensorProduct(x -> shapeFunction1
//				                          .gradientInCell(x, cell)
//				                          .inner(shapeFunction2.gradientInCell(x, cell))
//				                          * (Double) weight.value(x),
//			                          cell.getComponentCells(), QuadratureRule1D.Gauss5);
			return integrateOnReferenceCell(x -> shapeFunction1
				                                .gradientOnReferenceCell(x, cell)
				                                .inner(shapeFunction2.gradientOnReferenceCell(x, cell))
				                                * (Double) weight.value(x),
			                                cell);
		}
		throw new UnsupportedOperationException("unknown integral name");
	}
}
