package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface FastEvaluatedScalarShapeFunction<CT extends Cell<CT,FT>, FT extends Face<CT,FT>> extends FastEvaluatedScalarFunction, ScalarShapeFunction<CT,FT
	>
{
	double fastValueInCell(CoordinateVector pos, CT cell);
	double[] fastGradientInCell(CoordinateVector pos, CT cell);
	default double[][] fastHessianInCell(CoordinateVector pos, CT cell)
	{
		throw new UnsupportedOperationException();
	}
	
	@Override
	default Double valueInCell(CoordinateVector pos, CT cell)
	{
		return fastValueInCell(pos, cell);
	}
	
	@Override
	default CoordinateVector gradientInCell(CoordinateVector pos, CT cell)
	{
		return CoordinateVector.fromValues(fastGradientInCell(pos, cell));
	}
	
	@Override
	default CoordinateMatrix hessianInCell(CoordinateVector pos, CT cell)
	{
		return new CoordinateMatrix(fastHessianInCell(pos, cell));
	}
	
	@Override
	default double fastValue(CoordinateVector pos)
	{
		for(CT cell: getCells())
			if(cell.isInCell(pos))
				return fastValueInCell(pos, cell);
		return 0.;
	}
	
	@Override
	default double[] fastGradient(CoordinateVector pos)
	{
		for(CT cell: getCells())
			if(cell.isInCell(pos))
				return fastGradientInCell(pos, cell);
		return new double[getDomainDimension()];
	}
	
	@Override
	default double[][] fastHessian(CoordinateVector pos)
	{
		for(CT cell: getCells())
			if(cell.isInCell(pos))
				return fastHessianInCell(pos, cell);
		return new double[getDomainDimension()][getDomainDimension()];
	}
	
	@Override
	default Double value(CoordinateVector pos)
	{
		return fastValue(pos);
	}
	
	@Override
	default CoordinateVector gradient(CoordinateVector pos)
	{
		return new CoordinateVector(fastGradient(pos));
	}
	
	@Override
	default CoordinateMatrix hessian(CoordinateVector pos)
	{
		return new CoordinateMatrix(fastHessian(pos));
	}
}
