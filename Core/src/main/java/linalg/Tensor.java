package linalg;

import com.google.common.primitives.Ints;

import java.util.*;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public interface Tensor
{
	double at(int... coordinates);
	void set(double value, int... coordinates);
	void add(double value, int...coordinates);
	
	Tensor add(Tensor other);
	default Tensor sub(Tensor other)
	{
		if(Ints.toArray(getShape()) != Ints.toArray(other.getShape()))
			throw new IllegalArgumentException("Tensors are of different size");
		return add(other.mul(-1.));
	}
	Tensor mul(double scalar);
	
	Map<List<Integer>,Double> getCoordinateEntryList();
	List<? extends Tensor> unfoldDimension(int dimension);
	
	int getSparseEntryCount();
	boolean isSparse();
	
	int getOrder();
	List<Integer> getShape();
	default long size()
	{
		Optional<Long> ret = getShape().stream().map(Integer::longValue).reduce(Math::multiplyExact);
		if(ret.isPresent())
			return ret.get();
		else
			return 0;
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
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Tensors are of different size");
		if(tol.length > 1)
			throw new IllegalArgumentException("Only one Tolerance accepted");
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