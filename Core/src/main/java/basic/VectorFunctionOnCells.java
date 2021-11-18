package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

import java.util.function.BiFunction;

public interface VectorFunctionOnCells<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	extends FunctionOnCells<CT,
	FT, CoordinateVector, CoordinateMatrix, CoordinateTensor>, VectorFunction
{
	
	static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	VectorFunctionOnCells<CT, FT>
	concatenate(final VectorFunction outer, final VectorFunctionOnCells<CT, FT> inner)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (inner.getRangeDimension() != outer.getDomainDimension())
				throw new IllegalArgumentException("Inner function has the wrong range");
		}
		
		final VectorFunction function = outer.concatenateWith(inner);
		return new VectorFunctionOnCells<>()
		{
			@Override
			public CoordinateVector valueInCell(final CoordinateVector pos, final CT cell)
			{
				return function.value(pos);
			}
			
			@Override
			public CoordinateMatrix gradientInCell(final CoordinateVector pos, final CT cell)
			{
				return function.gradient(pos);
			}
			
			@Override
			public int getDomainDimension()
			{
				return function.getDomainDimension();
			}
			
			@Override
			public int getRangeDimension()
			{
				return function.getRangeDimension();
			}
			
			@Override
			public CoordinateVector defaultValue()
			{
				return function.defaultValue();
			}
			
			@Override
			public CoordinateMatrix defaultGradient()
			{
				return function.defaultGradient();
			}
			
			@Override
			public CoordinateTensor defaultHessian()
			{
				return function.defaultHessian();
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return function.value(pos);
			}
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return function.gradient(pos);
			}
			
			@Override
			public VectorFunction concatenateWith(final VectorFunction f)
			{
				return function.concatenateWith(f);
			}
		};
	}
	
	static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>> VectorFunctionOnCells<CT, FT>
	fromLambda(final java.util.function.Function<CoordinateVector, CoordinateVector> lambdaFunction,
	           final BiFunction<CoordinateVector, CT, CoordinateVector> lambdaFunctionOnCells,
	           final int domainDimension,
	           final int rangeDimension)
	{
		return new VectorFunctionOnCells<CT, FT>()
		{
			@Override
			public CoordinateVector valueInCell(final CoordinateVector pos, final CT cell)
			{
				return lambdaFunctionOnCells.apply(pos, cell);
			}
			
			@Override
			public int getRangeDimension()
			{
				return rangeDimension;
			}
			
			@Override
			public int getDomainDimension()
			{
				return domainDimension;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return lambdaFunction.apply(pos);
			}
		};
	}
}
