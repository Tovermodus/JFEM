package linalg;

import basic.PerformanceArguments;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Streams;
import com.google.common.primitives.Ints;

import java.util.*;
import java.util.stream.IntStream;

public interface Vector extends Tensor
{
	default int getLength()
	{
		return getShape().size();
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
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Vectors are of different size");
		return add(other.mul(-1.));
	}
	
	@Override
	Vector mul(double scalar);
	
	Matrix outer(Vector other);
	
	default double inner(Vector other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getLength() != other.getLength())
				throw new IllegalArgumentException("Vectors are of different size");
		return getShape().range().stream().mapToDouble(c -> at(c) * other.at(c)).sum();
	}
	
	default double euclidianNorm()
	{
		return Math.sqrt(this.inner(this));
	}
	
	default Vector vmMul(Matrix matrix)
	{
		return matrix.tvMul(this);
	}
	
	@Override
	default String printFormatted(double... tol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (tol.length > 1)
				throw new IllegalArgumentException("Only one Tolerance accepted");
		if (tol.length == 0)
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
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Vectors are of different size");
			if (tol.length > 1)
				throw new IllegalArgumentException("Only one Tolerance accepted");
		}
		double[] finalTol;
		if (tol.length == 0)
			finalTol = new double[]{0};
		else
			finalTol = tol;
		double absmax = absMaxElement() + other.absMaxElement();
		if (!(other instanceof Vector))
			return false;
		return getShape().range().stream().allMatch(c -> Math.abs(at(c) - other.at(c)) < finalTol[0] * absmax);
		
	}
	
}
