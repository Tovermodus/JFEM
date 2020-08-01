package linalg;

import com.google.common.primitives.Ints;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.IntStream;

public interface Vector extends Tensor
{
	@Override
	default int getOrder()
	{
		return 1;
	}
	@Override
	default Map<List<Integer>,Double> getCoordinateEntryList()
	{
		Map<List<Integer>,Double> ret = new TreeMap<>(new CoordinateComparator());
		for(int i = 0; i < getLength(); i++)
			if(at(i) != 0)
				ret.put(Ints.asList(i), at(i));
		return ret;
	}
	default int getLength()
	{
		return getShape().get(0);
	}
	@Override
	default List<Tensor> unfoldDimension(int dimension)
	{
		throw new UnsupportedOperationException("Vector can't unfold");
	}
	@Override
	Vector add(Tensor other);
	@Override
	default Vector sub(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		return add(other.mul(-1.));
	}
	
	@Override
	default void addInPlace(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		for (int i = 0; i < getLength(); i++)
		{
			set(at(i)+other.at(i),i);
		}
	}
	
	@Override
	default void subInPlace(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		for (int i = 0; i < getLength(); i++)
		{
			set(at(i)-other.at(i),i);
		}
	}
	
	@Override
	Vector mul(double scalar);
	
	default double inner(Vector other)
	{
		if(getLength() != other.getLength())
			throw new IllegalArgumentException("Vectors are of different size");
		return IntStream.range(0,getLength()).mapToDouble(i->at(i)*other.at(i)).sum();
	}
	default double euclidianNorm()
	{
		return Math.sqrt(this.inner(this));
	}
	default Vector vmMul(Matrix matrix)
	{
		return matrix.tvMul(this);
	}
//	{
//		assert getShape().get(0).equals(matrix.getShape().get(0));
//		return IntStream.range(0,matrix.getShape().get(1)).parallel().mapToDouble(i->matrix.unfoldDimension(1).get(i).inner(this)).sum();
//	}
	@Override
	default String printFormatted(double...tol)
	{
		if(tol.length > 1)
			throw new IllegalArgumentException("Only one Tolerance accepted");
		if(tol.length == 0)
			tol = new double[]{1e-14};
		String ret = "[";
		for(int i = 0; i < getLength(); i++)
		{
			if (Math.abs(at(i)) < tol[0])
				ret = ret.concat("<>........  ");
			else
			{
				if(at(i)>=0)
					ret = ret.concat("+");
				ret = ret.concat(String.format("%6.3e", at(i)) + "  ");
			}
		}
		ret = ret.concat("]");
		return ret;
	}
	@Override
	default boolean almostEqual(Tensor other, double... tol)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		if(tol.length > 1)
			throw new IllegalArgumentException("Only one Tolerance accepted");
		double[] finalTol;
		if(tol.length == 0)
			finalTol = new double[]{0};
		else
			finalTol=tol;
		double absmax = absMaxElement() + other.absMaxElement();
		if(!(other instanceof Vector))
			return false;
		return IntStream.range(0,getLength()).parallel().allMatch(i -> Math.abs(at(i) - other.at(i)) < finalTol[0] + absmax);
		
	}
	
}
