package basic;

import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public interface ScalarFunctionOnCells<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	extends FunctionOnCells<CT,
	FT, Double, CoordinateVector, CoordinateMatrix>, ScalarFunction
{
	
	static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>> ScalarFunctionOnCells<CT, FT> constantFunction(final double constant)
	{
		return new ScalarFunctionOnCells<>()
		{
			@Override
			public Double valueInCell(final CoordinateVector pos, final CT cell)
			{
				return constant;
			}
			
			@Override
			public CoordinateVector gradientInCell(final CoordinateVector pos, final CT cell)
			{
				return new CoordinateVector((int) pos.size());
			}
			
			@Override
			public int getDomainDimension()
			{
				return -1;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return constant;
			}
			
			@Override
			public CoordinateVector gradient(final CoordinateVector pos)
			{
				return new CoordinateVector((int) pos.size());
			}
			
			@Override
			public CoordinateMatrix hessian(final CoordinateVector pos)
			{
				return new CoordinateDenseMatrix((int) pos.size(), (int) pos.size());
			}
		};
	}
	
	static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	ScalarFunctionOnCells<CT, FT>
	concatenate(final ScalarFunction outer, final VectorFunctionOnCells<CT, FT> inner)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (inner.getRangeDimension() != outer.getDomainDimension())
				throw new IllegalArgumentException("Inner function has the wrong range");
		}
		
		final ScalarFunction function = outer.concatenateWith(inner);
		return new ScalarFunctionOnCells<CT, FT>()
		{
			@Override
			public Double valueInCell(final CoordinateVector pos, final CT cell)
			{
				return outer.value(inner.valueInCell(pos, cell));
			}
			
			@Override
			public CoordinateVector gradientInCell(final CoordinateVector pos, final CT cell)
			{
				return inner.gradientInCell(pos, cell)
				            .mvMul(outer.gradient(inner.valueInCell(pos, cell)));
			}
			
			@Override
			public int getDomainDimension()
			{
				return function.getDomainDimension();
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return function.value(pos);
			}
			
			@Override
			public CoordinateVector gradient(final CoordinateVector pos)
			{
				return function.gradient(pos);
			}
		};
	}
}
