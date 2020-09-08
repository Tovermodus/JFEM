package linalg;

import com.google.common.primitives.Ints;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public interface Matrix extends Tensor, VectorMultiplyable
{
	@Override
	default int getOrder()
	{
		return 2;
	}
	@Override
	Matrix add(Tensor other);
	@Override
	default Matrix sub(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Matrices are of different size");
		return add(other.mul(-1.));
	}
	
	@Override
	default void addInPlace(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Matrices are of different size");
		for(int i = 0; i < getRows(); i++)
			for (int j = 0; j < getCols(); j++)
			{
				set(at(i,j)+other.at(i,j),i,j);
			}
	}
	
	@Override
	default void subInPlace(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Matrices are of different size");
		for(int i = 0; i < getRows(); i++)
			for (int j = 0; j < getCols(); j++)
			{
				set(at(i,j)-other.at(i,j),i,j);
			}
	}
	
	@Override
	default Map<List<Integer>,Double> getCoordinateEntryList()
	{
		Map<List<Integer>,Double> ret = new TreeMap<>(new CoordinateComparator());
		for(int i = 0; i < getRows(); i++)
			for(int j = 0; j < getCols(); j++)
				if(at(i,j) != 0)
					ret.put(Ints.asList(i,j), at(i,j));
		return ret;
	}
	@Override
	Matrix mul(double scalar);
	
	default int getRows()
	{
		return getShape().get(0);
	}
	default int getCols()
	{
		return getShape().get(1);
	}
	@Override
	List<Vector> unfoldDimension(int dimension);
	
	Matrix transpose();
	
	@Override
	Vector mvMul(Vector vector);
	default Vector tvMul(Vector vector)
	{
		if(getRows()!=vector.getLength())
			throw new IllegalArgumentException("Incompatible sizes");
		return transpose().mvMul(vector);
	}
	default
	Double frobeniusInner(Matrix other)
	{
		throw new UnsupportedOperationException("not yet implemented");
	}
	Matrix mmMul(Matrix matrix);
	default Matrix tmMul(Matrix matrix)
	{
		if(getRows()!=(matrix.getRows()))
			throw new IllegalArgumentException("Incompatible sizes");
		return transpose().mmMul(matrix);
	}
	default Matrix mtMul(Matrix matrix)
	{
		if(getCols()!=(matrix.getCols()))
			throw new IllegalArgumentException("Incompatible sizes");
		return mmMul(matrix.transpose());
	}
	@Override
	default String printFormatted(double...tol)
	{
		if(tol.length > 1)
			throw new IllegalArgumentException("Only one Tolerance accepted");
		if(tol.length == 0)
			tol = new double[]{1e-14};
		String ret = "[";
		for(int i = 0; i < getRows(); i++)
		{
			for(int j = 0; j < getCols(); j++)
			{
				if (Math.abs(at(i,j)) < tol[0])
					ret = ret.concat("<>........  ");
				else
				{
					if (at(i,j) >= 0)
						ret = ret.concat("+");
					ret = ret.concat(String.format("%6.3e", at(i,j)) + "  ");
				}
			}
			ret = ret.concat("\n ");
		}
		ret = ret.concat("]");
		return ret;
	}
}
