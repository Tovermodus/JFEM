package linalg;

import com.google.common.primitives.Ints;

import java.util.*;

public class SparseVector implements Vector
{
	TreeMap<List<Integer>, Double> values;
	int size;
	public SparseVector(int size){
		this.size = size;
		values = new TreeMap<>(new CoordinateComparator());
	}
	public SparseVector(SparseVector vector){
		this.size = vector.size;
		values = new TreeMap<>(vector.values);
	}
	@Override
	public double at(int... coordinates)
	{
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if(values.containsKey(Ints.asList(coordinates)))
			return values.get(Ints.asList(coordinates));
		return 0;
	}
	
	@Override
	public synchronized void set(double value, int... coordinates)
	{
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		values.put(Ints.asList(coordinates), value);
		if(value == 0)
			values.remove(Ints.asList(coordinates));
		
	}
	
	@Override
	public synchronized void add(double value, int... coordinates)
	{
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		value += at(coordinates);
		values.put(Ints.asList(coordinates), value);
		if(value == 0)
			values.remove(Ints.asList(coordinates));
	}
	
	@Override
	public Vector add(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		if(other.isSparse())
		{
			SparseVector ret = new SparseVector(this);
			for(List<Integer> key: other.getCoordinateEntryList().keySet())
				ret.add(other.getCoordinateEntryList().get(key), Ints.toArray(key));
			return ret;
		}
		else
			return (Vector)other.add(this);
	}
	
	@Override
	public Vector sub(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		if(other.isSparse())
		{
			SparseVector ret = new SparseVector(this);
			for(List<Integer> key: other.getCoordinateEntryList().keySet())
				ret.add(-other.getCoordinateEntryList().get(key), Ints.toArray(key));
			return ret;
		}
		else
			return (Vector) other.sub(this).mul(-1);
	}
	@Override
	public Vector mul(double scalar)
	{
		SparseVector ret = new SparseVector(this.size);
		for(List<Integer> key: getCoordinateEntryList().keySet())
			ret.add(scalar* getCoordinateEntryList().get(key), Ints.toArray(key));
		return ret;
	}
	
	@Override
	public double inner(Vector other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		return values.keySet().stream().parallel().mapToDouble(key->values.get(key)*other.at(Ints.toArray(key))).sum();
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return values.size();
	}
	
	@Override
	public boolean isSparse()
	{
		return true;
	}
	
	@Override
	public List<Integer> getShape()
	{
return Ints.asList(size);
	}
	
	
	@Override
	public boolean equals(Object obj)
	{
		if(!(obj instanceof Vector))
			return false;
		return almostEqual((Tensor) obj);
	}
	
	@Override
	public String toString()
	{
		return printFormatted();
	}
//	@Override
//	public Vector vmMul(Matrix matrix)
//	{
//		if(!getShape().get(0).equals(matrix.getShape().get(0)))
//			throw new IllegalArgumentException("Incompatible sizes");
//		if(!matrix.isSparse())
//		{
//			DenseVector ret = new DenseVector(matrix.getShape().get(1));
//			ret.entries =
//				IntStream.range(0, matrix.getShape().get(1)).parallel().mapToDouble(i -> matrix.unfoldDimension(1).get(i).inner(this)).toArray();
//			return ret;
//		}
//	}
}
