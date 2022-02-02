package linalg;

import basic.PerformanceArguments;
import org.jetbrains.annotations.NotNull;

public class CoordinateVector
	extends DenseVector
	implements Comparable<CoordinateVector>
{
	private static final CoordinateVector unit21 = CoordinateVector.fromValues(1, 0);
	private static final CoordinateVector unit22 = CoordinateVector.fromValues(0, 1);
	private static final CoordinateVector zero = CoordinateVector.fromValues(0, 0);
	
	public CoordinateVector(final Vector v)
	{
		super(v);
		if (PerformanceArguments.getInstance().executeChecks)
			if (v.getLength() > 4) throw new IllegalArgumentException("only 1D, 2D, 3D and 4D supported");
	}
	
	public CoordinateVector(final int d)
	{
		super(d);
		if (PerformanceArguments.getInstance().executeChecks)
			if (d > 4) throw new IllegalArgumentException("only 1D, 2D, 3D and 4D supported");
	}
	
	public CoordinateVector(final double[] vector)
	{
		super(vector);
	}
	
	public static CoordinateVector getUnitVector(final int d, final int index)
	{
		if (d == 2)
		{
			if (index == 0)
				return unit21;
			else return unit22;
		}
		final CoordinateVector ret = new CoordinateVector(d);
		ret.set(1, index);
		return ret;
	}
	
	public double dist(final CoordinateVector other)
	{
		double sumOfSquares = 0;
		for (int i = 0; i < getLength(); i++)
			sumOfSquares += Math.pow(entries[i] - other.entries[i], 2);
		return Math.sqrt(sumOfSquares);
	}
	
	public static CoordinateVector getUnitVector(final int d, final int index, final double scale)
	{
		final CoordinateVector ret = new CoordinateVector(d);
		ret.set(scale, index);
		return ret;
	}
	
	public static CoordinateVector getZero(final int d)
	{
		if (d == 2)
			return zero;
		return new CoordinateVector(d);
	}
	
	public static CoordinateVector fromValues(final double... values)
	{
		final CoordinateVector ret = new CoordinateVector(values.length);
		for (int i = 0; i < values.length; i++)
		{
			ret.set(values[i], i);
		}
		return ret;
	}
	
	public CoordinateVector normalize()
	{
		return this.mul(1. / euclidianNorm());
	}
	
	public static CoordinateVector repeat(final double val, final int number)
	{
		final CoordinateVector ret = new CoordinateVector(number);
		for (int i = 0; i < number; i++)
		{
			ret.set(val, i);
		}
		return ret;
	}
	
	public CoordinateVector addCoordinate(final double time)
	{
		final CoordinateVector ret = new CoordinateVector(getLength() + 1);
		for (int i = 0; i < getLength(); i++)
		{
			ret.set(at(i), i);
		}
		ret.set(time, getLength());
		return ret;
	}
	
	public Rank1CoordinateMatrix outer(final CoordinateVector other)
	{
		return new Rank1CoordinateMatrix(this, other);
//		final CoordinateMatrix ret = new CoordinateMatrix(getLength(), other.getLength());
//		for (int i = 0; i < getLength(); i++)
//		{
//			for (int j = 0; j < other.getLength(); j++)
//			{
//				ret.set(at(i) * other.at(j), i, j);
//			}
//		}
//		return ret;
	}
	
	public CoordinateMatrix asCoordinateMatrix()
	{
		final CoordinateDenseMatrix ret = new CoordinateDenseMatrix(entries.length, 1);
		for (int i = 0; i < getLength(); i++)
		{
			ret.set(entries[i], i, 0);
		}
		return ret;
	}
	
	@Override
	public CoordinateVector add(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (getLength() != ((Vector) other).getLength())
				throw new IllegalArgumentException("Vectors are of different size");
		}
		final CoordinateVector ret = new CoordinateVector(this);
		if (!(other instanceof CoordinateVector))
		{
			for (int i = 0; i < getLength(); i++)
				ret.entries[i] = entries[i] + other.at(i);
		} else
		{
			for (int i = 0; i < getLength(); i++)
				ret.entries[i] = entries[i] + ((DenseVector) other).entries[i];
		}
		return ret;
	}
	
	@Override
	public CoordinateVector sub(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (getLength() != ((Vector) other).getLength())
				throw new IllegalArgumentException("Vectors are of different size");
		}
		final CoordinateVector ret = new CoordinateVector((Vector) other);
		for (int i = 0; i < getLength(); i++)
			ret.entries[i] = entries[i] - ret.entries[i];
		return ret;
	}
	
	@Override
	public CoordinateVector mul(final double scalar)
	{
		final CoordinateVector ret = new CoordinateVector(entries.length);
		for (int i = 0; i < getLength(); i++)
		{
			ret.set(at(i) * scalar, i);
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
		return entries[entries.length - 1];
	}
	
	@Override
	public int compareTo(@NotNull final CoordinateVector o)
	{
		return CoordinateComparator.comp(this, o);
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (obj instanceof CoordinateVector) return compareTo((CoordinateVector) obj) == 0;
		return false;
	}
	
	public CoordinateVector cross(final CoordinateVector other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (getLength() != 3 || other.getLength() != 3)
				throw new IllegalArgumentException("can't compute cross product");
		}
		return CoordinateVector.fromValues(entries[1] * other.entries[2] - entries[2] * other.entries[1],
		                                   entries[2] * other.entries[0] - entries[0] * other.entries[2],
		                                   entries[0] * other.entries[1] - entries[1] * other.entries[0]);
	}
}
