package distorted;

import basic.CellIntegral;
import basic.Function;
import basic.ScalarFunction;
import distorted.geometry.DistortedCell;
import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import tensorproduct.QuadratureRule1D;
import tensorproduct.geometry.TPCell;

import java.util.function.ToDoubleFunction;

public class DistortedVectorCellIntegral
	extends CellIntegral<DistortedCell, DistortedVectorShapeFunction>
{
	
	public static final String H1 = "H1";
	public static final String VALUE_VALUE = "ValueValue";
	public static final String GRAD_GRAD = "GradGrad";
	public static final String SYM_SYM = "SymGrad";
	public static final String TRACE_SYM = "TraceSym";
	
	public DistortedVectorCellIntegral(final double weight,
	                                   final String name,
	                                   final QuadratureRule1D quadratureRule1D)
	{
		this(ScalarFunction.constantFunction(weight), name, quadratureRule1D);
	}
	
	public DistortedVectorCellIntegral(final double weight, final String name)
	{
		this(ScalarFunction.constantFunction(weight), name);
	}
	
	public DistortedVectorCellIntegral(final Function<?, ?, ?> weight, final String name)
	{
		this(weight, name, QuadratureRule1D.Gauss5);
	}
	
	public DistortedVectorCellIntegral(final Function<?, ?, ?> weight,
	                                   final String name,
	                                   final QuadratureRule1D quadratureRule1D)
	{
		super(weight, name, quadratureRule1D);
		if (name.equals(H1) && !(weight.defaultValue() instanceof Double)) throw new IllegalArgumentException();
	}
	
	public DistortedVectorCellIntegral(final String name)
	{
		this(name, QuadratureRule1D.Gauss5);
	}
	
	public DistortedVectorCellIntegral(final String name, final QuadratureRule1D quadratureRule1D)
	{
		super(name, quadratureRule1D);
	}
	
	public static double integrateOnReferenceCell(final ToDoubleFunction<CoordinateVector> eval,
	                                              final DistortedCell cell,
	                                              final QuadratureRule1D quadratureRule)
	{
		double ret = 0;
		double val;
		final TPCell referenceCell = cell.referenceCell;
		final CoordinateVector quadraturePoint = new CoordinateVector(cell.getDimension());
		final double[][][] pointsWeights = new double[cell.getDimension()][2][quadratureRule.length()];
		for (int j = 0; j < cell.getDimension(); j++)
			pointsWeights[j] = referenceCell.getComponentCell(j)
			                                .distributeQuadrature(quadratureRule);
		final IntCoordinates quadraturePointSize = IntCoordinates.repeat(quadratureRule.length(),
		                                                                 cell.getDimension());
		for (final IntCoordinates c : quadraturePointSize.range())
		{
			for (int i = 0; i < c.getDimension(); i++)
			{
				quadraturePoint.set(pointsWeights[i][0][c.get(i)], i);
			}
			val = eval.applyAsDouble(quadraturePoint) * Math.abs(
				cell.transformationGradientFromReferenceCell(quadraturePoint)
				    .determinant());
			for (int i = 0; i < cell.getDimension(); i++)
				val *= pointsWeights[i][1][c.get(i)];
			ret += val;
		}
		return ret;
	}
	
	public static double integrateOnCell(final ToDoubleFunction<CoordinateVector> eval,
	                                     final DistortedCell cell,
	                                     final QuadratureRule1D quadratureRule)
	{
		double ret = 0;
		double val;
		final TPCell referenceCell = cell.referenceCell;
		final CoordinateVector quadraturePoint = new CoordinateVector(cell.getDimension());
		final double[][][] pointsWeights = new double[cell.getDimension()][2][quadratureRule.length()];
		for (int j = 0; j < cell.getDimension(); j++)
			pointsWeights[j] = referenceCell.getComponentCell(j)
			                                .distributeQuadrature(quadratureRule);
		final IntCoordinates quadraturePointSize = IntCoordinates.repeat(quadratureRule.length(),
		                                                                 cell.getDimension());
		for (final IntCoordinates c : quadraturePointSize.range())
		{
			for (int i = 0; i < c.getDimension(); i++)
			{
				quadraturePoint.set(pointsWeights[i][0][c.get(i)], i);
			}
			final CoordinateVector pointOnGrid = cell.transformFromReferenceCell(quadraturePoint);
			val = eval.applyAsDouble(pointOnGrid);
			for (int i = 0; i < cell.getDimension(); i++)
				val *= pointsWeights[i][1][c.get(i)];
			ret += val;
		}
		return ret;
	}
	
	@Override
	public double evaluateCellIntegral(final DistortedCell cell,
	                                   final DistortedVectorShapeFunction shapeFunction1,
	                                   final DistortedVectorShapeFunction shapeFunction2)
	{
		if (name.equals(SYM_SYM))
		{
			return integrateOnReferenceCell(x ->
			                                {
				                                final CoordinateMatrix grad1
					                                = shapeFunction1.gradientOnReferenceCell(x,
					                                                                         cell);
				                                final CoordinateMatrix grad2
					                                = shapeFunction2.gradientOnReferenceCell(x,
					                                                                         cell);
				
				                                return grad1
					                                .add(grad1.transpose())
					                                .frobeniusInner(grad2.add(
						                                grad2.transpose())) * (Double) weight.value(
					                                x) / 4;
			                                }, cell, quadratureRule1D);
		}
		if (name.equals(TRACE_SYM))
		{
			return integrateOnReferenceCell(x ->
			                                {
				                                final CoordinateMatrix grad1
					                                = shapeFunction1.gradientOnReferenceCell(x,
					                                                                         cell);
				                                final CoordinateMatrix grad2
					                                = shapeFunction2.gradientOnReferenceCell(x,
					                                                                         cell);
				                                final CoordinateMatrix traceId =
					                                CoordinateDenseMatrix.identity(x.getLength())
					                                                     .mul(grad1.add(grad1.transpose())
					                                                               .trace() / 4);
				                                return traceId
					                                .frobeniusInner(grad2.add(
						                                grad2.transpose())) * (Double) weight.value(
					                                x);
			                                }, cell, quadratureRule1D);
		}
		if (name.equals(H1))
		{
			return integrateOnReferenceCell(x -> (shapeFunction1
				                                      .valueOnReferenceCell(x, cell)
				                                      .inner(shapeFunction2.valueOnReferenceCell(x,
				                                                                                 cell)) + shapeFunction1
				                                .gradientOnReferenceCell(x, cell)
				                                .frobeniusInner(
					                                shapeFunction2.gradientOnReferenceCell(x, cell))) * (Double) weight.value(x),
			                                cell, quadratureRule1D);
		}
		if (name.equals(VALUE_VALUE))
		{
			return integrateOnReferenceCell(x -> (shapeFunction1
				                                      .valueOnReferenceCell(x, cell)
				                                      .inner(shapeFunction2.valueOnReferenceCell(x,
				                                                                                 cell)) * (Double) weight
				.value(x)), cell, quadratureRule1D);
		}
		if (name.equals(GRAD_GRAD))
		{
			return integrateOnReferenceCell(x -> (shapeFunction1
				                                      .gradientOnReferenceCell(x, cell)
				                                      .frobeniusInner(
					                                      shapeFunction2.gradientOnReferenceCell(x,
					                                                                             cell))) * (Double) weight
				.value(x), cell, quadratureRule1D);
		} else throw new IllegalArgumentException("Name unknown");
	}
}
