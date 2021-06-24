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
			getCoordinateEntryList().values().stream().parallel().mapToDouble(Double::doubleValue).max();
		if(max.isPresent())
			return max.getAsDouble();
		else
			return 0;
	}
	
	String printFormatted(double... tol);
	
	default boolean almostEqual(Tensor other, double... tol)
	{
		if(PerformanceArguments.getInstance().executeChecks)
		{
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Tensors are of different size");
			if (tol.length > 1)
				throw new IllegalArgumentException("Only one Tolerance accepted");
		}
		double[] finalTol;
		if(tol.length == 0)
			finalTol = new double[]{0};
		else
			finalTol=tol;
		double absmax = absMaxElement() + other.absMaxElement();
		List<? extends Tensor> unfolded = unfoldDimension(0);
		List<? extends Tensor> otherUnfolded = other.unfoldDimension(0);
		return IntStream.range(0,unfolded.size()).parallel().allMatch(i -> ((Tensor) unfolded.get(i)).almostEqual((Tensor) otherUnfolded.get(i),
			finalTol[0]+absmax*1e-14));
	}
}