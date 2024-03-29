package linalg;

import basic.DoubleCompare;
import basic.PerformanceArguments;
import org.apache.commons.math3.linear.RealVector;
import org.ujmp.core.Matrix;

import java.io.Serializable;
import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;

public class DenseVector
	implements MutableVector, MutableTensor, Serializable
{
	protected volatile double[] entries;
	
	public DenseVector(int size)
	{
		if (size <= 0) size = 1;
		entries = new double[size];
	}
	
	public DenseVector(final DenseVector vect, final boolean wrap)
	{
		if (wrap) entries = vect.entries;
		else entries = vect.entries.clone();
	}
	
	public DenseVector(final Vector vect)
	{
		entries = new double[vect.getLength()];
		for (int i = 0; i < entries.length; i++)
		{
			entries[i] = vect.at(i);
		}
	}
	
	public DenseVector(final double[] vect)
	{
		entries = vect.clone();
	}
	
	public DenseVector(final double[] vect, final boolean wrap)
	{
		if (wrap)
			entries = vect;
		else
			entries = vect.clone();
	}
	
	public static DenseVector vectorFromValues(final double... values)
	{
		final DenseVector ret = new DenseVector(values.length);
		for (int i = 0; i < values.length; i++)
			ret.set(values[i], i);
		return ret;
	}
	
	public static DenseVector getUnitVector(final int d, final int index)
	{
		final DenseVector ret = new DenseVector(d);
		ret.set(1, index);
		return ret;
	}
	
	public static DenseVector fromUJMPVector(final Matrix matrix)
	{
		final DenseVector ret = new DenseVector((int) (matrix.getRowCount() * matrix.getColumnCount()));
		for (int i = 0; i < ret.getLength(); i++)
			ret.set(matrix.getAsDouble(i, 0), i);
		return ret;
	}
	
	public static DenseVector fromRealVector(final RealVector vector)
	{
		final DenseVector ret = new DenseVector((int) (vector.getDimension()));
		for (int i = 0; i < ret.getLength(); i++)
			ret.set(vector.getEntry(i), i);
		return ret;
	}
	
	public static DenseVector concatenate(final Vector... vectors)
	{
		final DenseVector ret = new DenseVector(Arrays.stream(vectors)
		                                              .mapToInt(Vector::getLength)
		                                              .sum());
		int retIndex = 0;
		for (final Vector v : vectors)
			for (final IntCoordinates c : v.getShape()
			                               .range())
				ret.set(v.at(c), retIndex++);
		return ret;
	}
	
	public static DenseVector repeat(final double value, final int length)
	{
		final double[] values = new double[length];
		for (int i = 0; i < length; i++)
			values[i] = value;
		return DenseVector.vectorFromValues(values);
	}
	
	@Override
	public double inner(final Vector other)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (getLength() != other.getLength())
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
	public double absMaxElement()
	{
		double max = 0;
		for (int i = 0; i < getLength(); i++)
			if (Math.abs(entries[i]) > max) max = Math.abs(entries[i]);
		return max;
	}
	
	@Override
	public double at(final int... coordinates)
	{
		return entries[coordinates[0]];
	}
	
	@Override
	public void set(final double value, final int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 1) throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]] = value;
	}
	
	@Override
	public void add(final double value, final int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 1) throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]] += value;
	}
	
	@Override
	public void addInPlace(final Tensor other)
	{
		for (int i = 0; i < entries.length; i++)
			entries[i] += other.at(i);
	}
	
	@Override
	public void subInPlace(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (!other.getShape()
		                                                                .equals(getShape()))
			throw new IllegalArgumentException("Other has wrong shape");
		for (int i = 0; i < getLength(); i++)
			entries[i] -= other.at(i);
	}
	
	@Override
	public void mulInPlace(final double scalar)
	{
		for (int i = 0; i < entries.length; i++)
			entries[i] *= scalar;
	}
	
	@Override
	public int getLength()
	{
		return entries.length;
	}
	
	@Override
	public DenseVector add(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		final DenseVector ret = new DenseVector(this);
		if (!other.isSparse())
		{
			if (other instanceof DenseVector)
			{
				for (int i = 0; i < getLength(); i++)
					ret.entries[i] = entries[i] + ((DenseVector) other).entries[i];
			} else for (int i = 0; i < getLength(); i++)
			{
				ret.add(other.at(i), i);
			}
		} else for (final IntCoordinates key : other.getCoordinateEntryList()
		                                            .keySet())
			ret.add(other.getCoordinateEntryList()
			             .get(key), key);
		return ret;
	}
	
	@Override
	public double[] asArray()
	{
		return entries;
	}
	
	@Override
	public DenseVector sub(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		final DenseVector ret = new DenseVector(getLength());
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
	public DenseVector mul(final double scalar)
	{
		final DenseVector ret = new DenseVector(entries.length);
		for (int i = 0; i < getLength(); i++)
		{
			ret.entries[i] = entries[i] * scalar;
		}
		return ret;
	}
	
	public DenseVector componentWise(final DoubleUnaryOperator action)
	{
		final DenseVector ret = new DenseVector(this);
		for (int i = 0; i < getLength(); i++)
			ret.entries[i] = action.applyAsDouble(ret.entries[i]);
		return ret;
	}
	
	@Override
	public DenseMatrix outer(final Vector other)
	{
		final DenseMatrix ret = new DenseMatrix(getLength(), other.getLength());
		for (int i = 0; i < getLength(); i++)
		{
			for (int j = 0; j < other.getLength(); j++)
			{
				ret.set(at(i) * other.at(j), i, j);
			}
		}
		return ret;
	}
	
	public void addSmallVectorAt(final Vector small, final int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (coordinates.length != 1) throw new IllegalArgumentException("Wrong number of coordinates");
			if (coordinates[0] + small.getLength() > getLength())
				throw new IllegalArgumentException("small Vector too large position");
		}
		for (int i = 0; i < small.getLength(); i++)
		{
			add(small.at(i), i + coordinates[0]);
		}
	}
	
	no.uib.cipr.matrix.DenseVector toMTJvector()
	{
		final no.uib.cipr.matrix.DenseVector m = new no.uib.cipr.matrix.DenseVector(getShape().get(0));
		for (int i = 0; i < getLength(); i++)
		{
			m.set(i, at(i));
		}
		return m;
	}
	
	org.ujmp.core.Matrix toUJMPVector()
	{
		final org.ujmp.core.Matrix mat = org.ujmp.core.DenseMatrix.Factory.zeros(getLength(), 1);
		for (int i = 0; i < getLength(); i++)
			mat.setAsDouble(at(i), i, 0);
		return mat;
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (!(obj instanceof Vector)) return false;
		return almostEqual((Tensor) obj);
	}
	
	@Override
	public int hashCode()
	{
		int ret = 0;
		for (int i = 0; i < getLength(); i++)
		{
			final int elementHash = DoubleCompare.doubleHash(at(i));
			ret = ret * 37 + (i + 2) * (elementHash < 0 ? 3 : 2 + i + 3) * elementHash;
		}
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
	
	public double[] toArray()
	{
		return entries.clone();
	}
}
