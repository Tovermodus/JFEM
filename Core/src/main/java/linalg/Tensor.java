package linalg;

import basic.DoubleCompare;
import basic.PerformanceArguments;
import com.google.common.collect.ImmutableMap;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;

public interface Tensor
{
	double at(int... coordinates);
	
	default double at(final IntCoordinates coordinates)
	{
		return at(coordinates.asArray());
	}
	
	Tensor add(Tensor other);
	
	default Tensor sub(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getShape() != other.getShape())
				throw new IllegalArgumentException("Tensors are of different size");
		return add(other.mul(-1.));
	}
	
	Tensor mul(double scalar);
	
	default int getOrder()
	{
		return getShape().getDimension();
	}
	
	default ImmutableMap<IntCoordinates, Double> getCoordinateEntryList()
	{
		final Map<IntCoordinates, Double> ret = new HashMap<>();
		for (final IntCoordinates c : getShape().range())
			if (!DoubleCompare.almostEqual(at(c), 0))
				ret.put(c, at(c));
		return ImmutableMap.copyOf(ret);
	}
	
	List<? extends Tensor> unfoldDimension(int dimension);
	
	Tensor slice(IntCoordinates start, IntCoordinates end);
	
	int getSparseEntryCount();
	
	boolean isSparse();
	
	IntCoordinates getShape();
	
	default long size()
	{
		return getShape().size();
	}
	
	default double absMaxElement()
	{
		final OptionalDouble max =
			getCoordinateEntryList().values()
			                        .stream()
			                        .mapToDouble(Double::doubleValue)
			                        .map(
				                        Math::abs)
			                        .max();
		if (max.isPresent())
			return max.getAsDouble();
		else
			return 0;
	}
	
	String printFormatted(double... tol);
	
	default boolean almostZero()
	{
		return DoubleCompare.almostEqualAfterOps(absMaxElement(), 0, 0, size());
	}
	
	default boolean almostEqual(final Tensor other)
	{
		return almostEqual(other, 0);
	}
	
	default boolean almostEqual(final Tensor other, final double extraTol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Tensors are of different size" + getShape() + " " +
					                                   "!= " + other.getShape());
		}
		final double absmax = absMaxElement() + other.absMaxElement();
		for (final IntCoordinates c : getShape().range())
		{
			if (!DoubleCompare.almostEqualAfterOps(at(c),
			                                       other.at(c), absmax, size(), extraTol))
			{
				System.out.println(at(c) + " != " + other.at(c) +
					                   " difference: " + Math.abs(at(c) - other.at(c)) + " ");
				return false;
			}
		}
		return true;
	}
	
	default boolean almostEqualMute(final Tensor other)
	{
		return almostEqualMute(other, 1e-8);
	}
	
	default boolean almostEqualMute(final Tensor other, final double extraTol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Tensors are of different size");
		}
		final double absmax = absMaxElement() + other.absMaxElement();
		for (final IntCoordinates c : getShape().range())
		{
			if (!DoubleCompare.almostEqualAfterOps(at(c),
			                                       other.at(c), absmax, size(), extraTol))
			{
				return false;
			}
		}
		return true;
	}
}
