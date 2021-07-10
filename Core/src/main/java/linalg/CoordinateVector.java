package linalg;

import basic.PerformanceArguments;
import org.jetbrains.annotations.NotNull;

public class CoordinateVector extends DenseVector implements Comparable<CoordinateVector>
{
	public CoordinateVector(Vector v)
	{
		super(v);
		if(PerformanceArguments.getInstance().executeChecks)
		if(v.getLength()>4)
			throw new IllegalArgumentException("only 1D, 2D, 3D and 4D supported");
	}
	public CoordinateVector(int d)
	{
		super(d);
		if(PerformanceArguments.getInstance().executeChecks)
		if(d > 4)
			throw new IllegalArgumentException("only 1D, 2D, 3D and 4D supported");
	}
	public CoordinateVector(double [] vector)
	{
		super(vector);
	}
	public static CoordinateVector getUnitVector(int d, int index)
	{
		CoordinateVector ret = new CoordinateVector(d);
		ret.set(1,index);
		return ret;
	}
	public static CoordinateVector fromValues(double ... values)
	{
		CoordinateVector ret = new CoordinateVector(values.length);
		for (int i = 0; i < values.length; i++)
		{
			ret.set(values[i],i);
		}
		return ret;
	}
	
	public CoordinateVector addTime(double time)
	{
		CoordinateVector ret = new CoordinateVector(getLength()+1);
		for (int i = 0; i < getLength(); i++)
		{
			ret.set(at(i),i);
		}
		ret.set(time, getLength());
		return ret;
	}
	public CoordinateMatrix outer(CoordinateVector other)
	{
		return new CoordinateMatrix(super.outer(other));
	}
	public CoordinateMatrix asCoordinateMatrix()
	{
		CoordinateMatrix ret = new CoordinateMatrix(entries.length, 1);
		for (int i = 0; i < getLength(); i++)
		{
			ret.set(entries[i], i,0);
		}
		return ret;
	}
	@Override
	public CoordinateVector add(Tensor other)
	{
		if(PerformanceArguments.getInstance().executeChecks)
		{
			if (!(other instanceof CoordinateVector))
				throw new IllegalArgumentException("can't add coordinates and other vectors");
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Vectors are of different size");
		}
		CoordinateVector ret = new CoordinateVector(this);
		for (int i = 0; i < getLength(); i++)
			ret.entries[i] = entries[i]+((DenseVector) other).entries[i];
		return ret;
	}
	@Override
	public CoordinateVector sub(Tensor other)
	{
		if(PerformanceArguments.getInstance().executeChecks)
		{
			if (!(other instanceof CoordinateVector))
				throw new IllegalArgumentException("can't add coordinates and other vectors");
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Vectors are of different size");
		}
		CoordinateVector ret = new CoordinateVector(this);
		for (int i = 0; i < getLength(); i++)
			ret.entries[i] = entries[i]-((DenseVector) other).entries[i];
		return ret;
	}
	
	@Override
	public CoordinateVector mul(double scalar)
	{
		CoordinateVector ret = new CoordinateVector(entries.length);
		for (int i = 0; i < getLength(); i++)
		{
			ret.set(at(i)*scalar,i);
		}
		return ret;
	}
	
	public double x()
	{
		return entries[0];
	}
	public double y()
	{
		return entries[1];
	}
	public double z()
	{
		return entries[2];
	}
	public double t()
	{
		return entries[entries.length-1];
	}
	
	@Override
	public int compareTo(@NotNull CoordinateVector o)
	{
		return CoordinateComparator.comp(this, o);
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(obj instanceof CoordinateVector)
			return compareTo((CoordinateVector)obj) == 0;
		return false;
	}
}
