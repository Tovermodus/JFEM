package linalg;

import basic.PerformanceArguments;
import com.google.common.primitives.Ints;
import org.ujmp.core.Matrix;
import tensorproduct.TPEdge;

import java.util.List;

public class DenseVector implements MutableVector, MutableTensor
{
	protected volatile double[] entries;
	
	public DenseVector(int size)
	{
		if (size <= 0)
			size = 1;
		entries = new double[size];
	}
	
	public DenseVector(Vector vect)
	{
		entries = new double[vect.getLength()];
		for (int i = 0; i < vect.getLength(); i++)
		{
			entries[i] = vect.at(i);
		}
	}
	
	public DenseVector(double[] vect)
	{
		entries = vect.clone();
	}
	
	public static DenseVector vectorFromValues(double... values)
	{
		DenseVector ret = new DenseVector(values.length);
		for (int i = 0; i < values.length; i++)
			ret.set(values[i], i);
		return ret;
	}
	
	@Override
	public double inner(Vector other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getLength() != other.getLength())
				throw new IllegalArgumentException("Vectors have different size");
		if (other instanceof DenseVector)
		{
			double ret = 0;
			for (int i = 0; i < getLength(); i++)
				ret += this.entries[i] * ((DenseVector) other).entries[i];
			return ret;
		} else
		{
			double ret = 0;
			for (int i = 0; i < getLength(); i++)
				ret += this.entries[i] * other.at(i);
			return ret;
		}
		
		
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
	public IntCoordinates getShape()
	{
		return new IntCoordinates(entries.length);
	}
	
	@Override
	public double at(int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 1)
				throw new IllegalArgumentException("Wrong number of coordinates");
		return entries[coordinates[0]];
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 1)
				throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]] = value;
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 1)
				throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]] += value;
	}
	
	@Override
	public void addInPlace(Tensor other)
	{
		entries = add(other).entries;
	}
	
	@Override
	public void subInPlace(Tensor other)
	{
		entries = sub(other).entries;
	}
	
	@Override
	public void mulInPlace(double scalar)
	{
		entries = mul(scalar).entries;
	}
	
	@Override
	public DenseVector add(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Vectors are of different size");
		DenseVector ret = new DenseVector(this);
		if (!other.isSparse())
		{
			if (other instanceof DenseVector)
			{
				for (int i = 0; i < getLength(); i++)
					ret.entries[i] = entries[i] + ((DenseVector) other).entries[i];
			} else
				for (int i = 0; i < getLength(); i++)
				{
					ret.add(other.at(i), i);
				}
		} else
			for (IntCoordinates key : other.getCoordinateEntryList().keySet())
				ret.add(other.getCoordinateEntryList().get(key), key);
		return ret;
	}
	
	@Override
	public DenseVector sub(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Vectors are of different size");
		DenseVector ret = new DenseVector(getLength());
		if (other instanceof DenseVector)
		{
			for (int i = 0; i < getLength(); i++)
				ret.entries[i] = entries[i] - ((DenseVector) other).entries[i];
		} else
		{
			for (int i = 0; i < getLength(); i++)
			{
				ret.add(-other.at(i), i);
			}
		}
		return ret;
	}
	
	@Override
	public DenseVector mul(double scalar)
	{
		DenseVector ret = new DenseVector(entries.length);
		for (int i = 0; i < getLength(); i++)
		{
			ret.entries[i] = entries[i]*scalar;
		}
		return ret;
	}
	
	@Override
	public DenseMatrix outer(Vector other)
	{
		DenseMatrix ret = new DenseMatrix(getLength(), other.getLength());
		for (int i = 0; i < getLength(); i++)
		{
			for (int j = 0; j < other.getLength(); j++)
			{
				ret.set(at(i) * other.at(j), i, j);
			}
		}
		return ret;
	}
	
	no.uib.cipr.matrix.DenseVector toMTJvector()
	{
		no.uib.cipr.matrix.DenseVector m = new no.uib.cipr.matrix.DenseVector(getShape().get(0));
		for (int i = 0; i < getLength(); i++)
		{
			m.set(i, at(i));
		}
		return m;
	}
	
	org.ujmp.core.Matrix toUJMPVector()
	{
		org.ujmp.core.Matrix mat = org.ujmp.core.DenseMatrix.Factory.zeros(getLength(), 1);
		for (int i = 0; i < getLength(); i++)
			mat.setAsDouble(at(i), i, 0);
		return mat;
	}
	
	static DenseVector fromUJMPVector(Matrix matrix)
	{
		DenseVector ret = new DenseVector((int) (matrix.getRowCount() * matrix.getColumnCount()));
		for (int i = 0; i < ret.getLength(); i++)
			ret.set(matrix.getAsDouble(i, 0), i);
		return ret;
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if (!(obj instanceof Vector))
			return false;
		return almostEqual((Tensor) obj);
	}
	
	@Override
	public int hashCode()
	{
		int ret = -1024*entries.length;
		for(int i = 0; i < entries.length; i++)
			ret += entries[i]*Math.pow(2, i%30);
		return ret;
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
