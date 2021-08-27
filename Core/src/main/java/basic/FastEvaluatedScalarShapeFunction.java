package basic;

import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface FastEvaluatedScalarShapeFunction<CT extends Cell<CT, FT>, FT extends Face<CT, FT>> extends FastEvaluatedScalarFunction, ScalarShapeFunction<CT, FT>
{
	double fastValueInCell(CoordinateVector pos, CT cell);
	
	double[] fastGradientInCell(CoordinateVector pos, CT cell);
	
	default double[][] fastHessianInCell(final CoordinateVector pos, final CT cell)
	{
		throw new UnsupportedOperationException();
	}
	
	@Override
	default Double valueInCell(final CoordinateVector pos, final CT cell)
	{
		return fastValueInCell(pos, cell);
	}
	
	@Override
	default CoordinateVector gradientInCell(final CoordinateVector pos, final CT cell)
	{
		return CoordinateVector.fromValues(fastGradientInCell(pos, cell));
	}
	
	@Override
	default CoordinateMatrix hessianInCell(final CoordinateVector pos, final CT cell)
	{
		return new CoordinateDenseMatrix(fastHessianInCell(pos, cell));
	}
	
	@Override
	default double fastValue(final CoordinateVector pos)
	{
		for (final CT cell : getCells())
			if (cell.isInCell(pos)) return fastValueInCell(pos, cell);
		return 0.;
	}
	
	@Override
	default double[] fastGradient(final CoordinateVector pos)
	{
		for (final CT cell : getCells())
			if (cell.isInCell(pos)) return fastGradientInCell(pos, cell);
		return new double[getDomainDimension()];
	}
	
	@Override
	default double[][] fastHessian(final CoordinateVector pos)
	{
		for (final CT cell : getCells())
			if (cell.isInCell(pos)) return fastHessianInCell(pos, cell);
		return new double[getDomainDimension()][getDomainDimension()];
	}
	
	@Override
	default Double value(final CoordinateVector pos)
	{
		return fastValue(pos);
	}
	
	@Override
	default CoordinateVector gradient(final CoordinateVector pos)
	{
		return new CoordinateVector(fastGradient(pos));
	}
	
	@Override
	default CoordinateMatrix hessian(final CoordinateVector pos)
	{
		return new CoordinateDenseMatrix(fastHessian(pos));
	}
}
