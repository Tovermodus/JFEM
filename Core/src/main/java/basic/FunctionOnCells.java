package basic;

import linalg.CoordinateVector;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Stream;

public interface FunctionOnCells<CT extends Cell<CT, FT>, FT extends Face<CT, FT>, valueT, gradientT, hessianT>
	extends Function<valueT, gradientT, hessianT>
{
	valueT valueInCell(CoordinateVector pos, CT cell);
	
	default gradientT gradientInCell(final CoordinateVector pos, final CT cell)
	{
		throw new UnsupportedOperationException();
	}
	
	default hessianT hessianInCell(final CoordinateVector pos, final CT cell)
	{
		throw new UnsupportedOperationException();
	}
	
	default Map<CoordinateVector, valueT> valuesInPointsOnCells(final List<CoordinateVector> points,
	                                                            final Collection<CT> cells)
	{
		
		final ConcurrentHashMap<CoordinateVector, CT> pointToCell = new ConcurrentHashMap<>();
		points.stream()
		      .parallel()
		      .forEach(p -> cells.stream()
		                         .filter(c -> c.isInCell(p))
		                         .findAny()
		                         .ifPresent(c -> pointToCell.put(p, c)));
		final ConcurrentHashMap<CoordinateVector, valueT> ret = new ConcurrentHashMap<>();
		Stream<CoordinateVector> stream = points.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads) stream = stream.parallel();
		stream.forEach(point -> ret.put(point, valueInCell(point, pointToCell.get(point))));
		return ret;
	}
	
	static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>, valueT, gradientT, hessianT>
	FunctionOnCells<CT, FT, valueT, gradientT, hessianT>
	concatenate(final Function<valueT, gradientT, hessianT> outer, final VectorFunctionOnCells<CT, FT> inner)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (inner.getRangeDimension() != outer.getDomainDimension())
				throw new IllegalArgumentException("Inner function has the wrong range");
		}
		
		final Function<valueT, gradientT, hessianT> function = outer.concatenateWith(inner);
		return new FunctionOnCells<>()
		{
			@Override
			public valueT valueInCell(final CoordinateVector pos, final CT cell)
			{
				return function.value(pos);
			}
			
			@Override
			public gradientT gradientInCell(final CoordinateVector pos, final CT cell)
			{
				return function.gradient(pos);
			}
			
			@Override
			public int getDomainDimension()
			{
				return function.getDomainDimension();
			}
			
			@Override
			public valueT defaultValue()
			{
				return function.defaultValue();
			}
			
			@Override
			public gradientT defaultGradient()
			{
				return function.defaultGradient();
			}
			
			@Override
			public hessianT defaultHessian()
			{
				return function.defaultHessian();
			}
			
			@Override
			public valueT value(final CoordinateVector pos)
			{
				return function.value(pos);
			}
			
			@Override
			public gradientT gradient(final CoordinateVector pos)
			{
				return function.gradient(pos);
			}
			
			@Override
			public Function<valueT, gradientT, hessianT> concatenateWith(final VectorFunction f)
			{
				return function.concatenateWith(f);
			}
		};
	}
}
