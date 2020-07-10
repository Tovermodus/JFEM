package linalg;

import com.google.common.primitives.Ints;
import org.ujmp.core.Matrix;

import java.util.List;

public class DenseVector implements Vector
{
	protected double [] entries;
	public DenseVector(int size)
	{
		if(size <= 0)
			size = 1;
		entries = new double[size];
	}
	public DenseVector(Vector vect)
	{
		entries = new double[vect.getLength()];
		for(List<Integer> key: vect.getCoordinateEntryList().keySet())
			add(vect.getCoordinateEntryList().get(key), Ints.toArray(key));
	}
	
	public static DenseVector vectorFromValues(double... values)
	{
		DenseVector ret = new DenseVector(values.length);
		for(int i = 0; i < values.length; i++)
			ret.set(values[i],i);
		return ret;
	}
	@Override
	public int getSparseEntryCount()
	{
		return -1;
	}
	
	@Override
	public boolean isSparse()
	{
		return false;
	}
	
	@Override
	public List<Integer> getShape()
	{
		return Ints.asList(entries.length);
	}
	
	@Override
	public double at(int... coordinates)
	{
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		return entries[coordinates[0]];
	}
	
	@Override
	public synchronized void set(double value, int... coordinates)
	{
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]] = value;
	}
	
	@Override
	public synchronized void add(double value, int... coordinates)
	{
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]] += value;
	}
	
	@Override
	public Vector add(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		DenseVector ret = new DenseVector(this);
		if(!other.isSparse())
			for (int i = 0; i < getLength(); i++)
			{
				ret.add(other.at(i),i);
			}
		else
			for(List<Integer> key: other.getCoordinateEntryList().keySet())
				ret.add(other.getCoordinateEntryList().get(key), Ints.toArray(key));
		return ret;
	}
	@Override
	public Vector sub(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		DenseVector ret = new DenseVector(this);
		if(!other.isSparse())
			for (int i = 0; i < getLength(); i++)
			{
				ret.add(-other.at(i),i);
			}
		else
			for(List<Integer> key: other.getCoordinateEntryList().keySet())
				ret.add(-other.getCoordinateEntryList().get(key), Ints.toArray(key));
		return ret;
	}
	@Override
	public Vector mul(double scalar)
	{
		DenseVector ret = new DenseVector(entries.length);
		for (int i = 0; i < getLength(); i++)
		{
			ret.set(at(i)*scalar,i);
		}
		return ret;
	}
	no.uib.cipr.matrix.DenseVector toMTJvector()
	{
		no.uib.cipr.matrix.DenseVector m = new no.uib.cipr.matrix.DenseVector(getShape().get(0));
		for(int i = 0; i < getLength(); i++)
		{
			m.set(i,at(i));
		}
		return m;
	}
	org.ujmp.core.Matrix toUJMPVector()
	{
		org.ujmp.core.Matrix mat = org.ujmp.core.DenseMatrix.Factory.zeros(getLength(),1);
		for(int i = 0; i < getLength(); i++)
			mat.setAsDouble(at(i), i,0);
		return mat;
	}
	static DenseVector fromUJMPVector(Matrix matrix)
	{
		DenseVector ret = new DenseVector((int) (matrix.getRowCount()*matrix.getColumnCount()));
		for(int i = 0; i < ret.getLength(); i++)
			ret.set(matrix.getAsDouble(i,0),i);
		return ret;
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
	
	public double[] getEntries()
	{
		return entries;
	}
}
