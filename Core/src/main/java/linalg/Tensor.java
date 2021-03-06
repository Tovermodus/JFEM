package linalg;

import basic.PerformanceArguments;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.primitives.Ints;
import java.util.*;
import java.util.stream.IntStream;

public interface Tensor
{
	double at(int ... coordinates);
	default double at(IntCoordinates coordinates)
	{
		return at(coordinates.asArray());
	}
	
	Tensor add(Tensor other);
	default Tensor sub(Tensor other)
	{
		if(PerformanceArguments.getInstance().executeChecks)
			if (getShape() != other.getShape())
				throw new IllegalArgumentException("Tensors are of different size");
		return add(other.mul(-1.));
	}
	Tensor mul(double scalar);
	
	default int getOrder()
	{
		return getShape().getDimension();
	}
	default ImmutableMap<IntCoordinates,Double> getCoordinateEntryList(){
		Map<IntCoordinates,Double> ret = new HashMap<>();
		for(IntCoordinates c: getShape().range())
			if(at(c) != 0)
				ret.put(c, at(c));
		return ImmutableMap.copyOf(ret);
	}
	List<? extends Tensor> unfoldDimension(int dimension);
	
	int getSparseEntryCount();
	boolean isSparse();
	
	IntCoordinates getShape();
	default long size()
	{
		return getShape().size();
	}
	
	default double absMaxElement()
	{
		OptionalDouble max =
			getCoordinateEntryList().values().stream().parallel().mapToDouble(Double::doubleValue).map(Math::abs).max();
		if(max.isPresent())
			return max.getAsDouble();
		else
			return 0;
	}
	
	String printFormatted(double... tol);
	
	default boolean almostEqual(Tensor other)
	{
		return almostEqual(other, 1e-15);
	}
	default boolean almostEqual(Tensor other, double tol)
	{
		if(PerformanceArguments.getInstance().executeChecks)
		{
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Tensors are of different size");
		}
		double absmax = absMaxElement() + other.absMaxElement();
		for(IntCoordinates c: getShape().range())
		{
			if(Math.abs(at(c) - other.at(c)) > Math.abs(tol*(1+absmax)))
			{
				System.out.println(at(c)+" != " +other.at(c)+" with tol " +  tol*(1+absmax)+ " " +
					" difference: " + Math.abs(at(c) - other.at(c)) + " ");
				return false;
			}
		}
		return true;
	}
}