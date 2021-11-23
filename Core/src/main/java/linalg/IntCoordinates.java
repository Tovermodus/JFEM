package linalg;

import com.google.common.collect.Streams;
import org.jetbrains.annotations.NotNull;

import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;

public class IntCoordinates
	implements Comparable<IntCoordinates>, Cloneable
{
	final private int[] coordinates;
	
	public IntCoordinates(final int... coordinates)
	{
		this.coordinates = coordinates;
	}
	
	public IntCoordinates(final List<Integer> coordinates)
	{
		this.coordinates = coordinates.stream()
		                              .mapToInt(Integer::intValue)
		                              .toArray();
	}
	
	public IntCoordinates(final IntCoordinates source)
	{
		this.coordinates = source.coordinates.clone();
	}
	
	public static IntCoordinates repeat(final int number, final int times)
	{
		final int[] coords = new int[times];
		for (int i = 0; i < times; i++)
		{
			coords[i] = number;
		}
		return new IntCoordinates(coords);
	}
	
	public int[] asArray()
	{
		return coordinates.clone();
	}
	
	public int get(final int d)
	{
		return coordinates[d];
	}
	
	public int getDimension()
	{
		return coordinates.length;
	}
	
	public IntCoordinates concatenate(final int[] outer)
	{
		final int[] concat = new int[this.getDimension()];
		for (int i = 0; i < concat.length; i++)
		{
			if (get(i) < 0 || get(i) > outer.length)
				throw new IllegalArgumentException("inner index <0 or >length");
			concat[i] = outer[get(i)];
		}
		return new IntCoordinates(concat);
	}
	
	public IntCoordinates add(final IntCoordinates o)
	{
		final IntCoordinates ret = new IntCoordinates(this);
		for (int i = 0; i < getDimension(); i++)
			ret.coordinates[i] += o.coordinates[i];
		return ret;
	}
	
	public IntCoordinates get_modified(final int d, final int value)
	{
		final IntCoordinates ret = new IntCoordinates(this);
		ret.coordinates[d] = value;
		return ret;
	}
	
	public IntCoordinates sub(final IntCoordinates o)
	{
		final IntCoordinates ret = new IntCoordinates(this);
		for (int i = 0; i < getDimension(); i++)
			ret.coordinates[i] -= o.coordinates[i];
		return ret;
	}
	
	public IntCoordinates zerosLike()
	{
		final int[] zeros = new int[coordinates.length];
		for (int i = 0; i < getDimension(); i++)
			zeros[i] = 0;
		return new IntCoordinates(zeros);
	}
	
	private IntCoordinates increaseLastDimension()
	{
		final IntCoordinates ret = new IntCoordinates(this);
		ret.coordinates[getDimension() - 1]++;
		return ret;
	}
	
	private IntCoordinates wrap(final IntCoordinates lowerBounds, final IntCoordinates upperBounds)
	{
		final IntCoordinates ret = new IntCoordinates(this);
		for (int d = getDimension() - 1; d > 0; d--)
			if (ret.coordinates[d] == upperBounds.coordinates[d])
			{
				ret.coordinates[d] = lowerBounds.coordinates[d];
				ret.coordinates[d - 1]++;
			}
		return ret;
	}
	
	@Override
	public int hashCode()
	{
		int ret = 0;
		for (final int c : coordinates)
		{
			ret = ret % 1657877;
			ret *= 104729;
			ret += c;
		}
		return ret;
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (obj instanceof IntCoordinates)
		{
			if (((IntCoordinates) obj).getDimension() != getDimension())
				return false;
			return this.compareTo((IntCoordinates) obj) == 0;
		}
		return false;
	}
	
	@Override
	public String toString()
	{
		final StringBuilder ret = new StringBuilder();
		for (final int c : coordinates)
			ret.append(c)
			   .append(" ");
		return ret.toString();
	}
	
	@Override
	public int compareTo(@NotNull final IntCoordinates o)
	{
		
		if (this.coordinates.length != o.coordinates.length)
			throw new IllegalArgumentException("coordinates have different length");
		for (int i = 0; i < coordinates.length; i++)
		{
			if (coordinates[i] > o.coordinates[i])
				return 1;
			if (coordinates[i] < o.coordinates[i])
				return -1;
		}
		return 0;
	}
	
	public int size()
	{
		int ret = 1;
		for (int i = 0; i < getDimension(); i++)
			ret *= coordinates[i];
		return ret;
	}
	
	public Range range()
	{
		return new Range(this);
	}
	
	public static class Range
		implements Iterable<IntCoordinates>, Iterator<IntCoordinates>
	{
		final IntCoordinates lowerBounds;
		final IntCoordinates upperBounds;
		IntCoordinates pointer;
		
		private int distanceFromUpper()
		{
			int ret = 0;
			for (int i = 0; i < upperBounds.getDimension(); i++)
				ret += upperBounds.get(i) - pointer.get(i);
			return ret;
		}
		
		public Range(final int[] lowerBounds, final int[] upperBounds)
		{
			this(new IntCoordinates(lowerBounds), new IntCoordinates(upperBounds));
		}
		
		public Range(final IntCoordinates lowerBounds, final IntCoordinates upperBounds)
		{
			if (lowerBounds.getDimension() != upperBounds.getDimension())
				throw new IllegalArgumentException("Range: bounds must have same length!" + lowerBounds + " , " + upperBounds);
			for (int i = 0; i < lowerBounds.getDimension(); i++)
			{
				if (lowerBounds.asArray()[i] >= upperBounds.asArray()[i])
					throw new IllegalArgumentException("Range: lower index higher than upper! " +
						                                   "But difference must be at least one" + lowerBounds + " >= " + upperBounds);
			}
			this.lowerBounds = lowerBounds;
			this.upperBounds = upperBounds;
		}
		
		public Range(final IntCoordinates upperBounds)
		{
			this(upperBounds.zerosLike(), upperBounds);
		}
		
		public Stream<IntCoordinates> stream()
		{
			return Streams.stream((Iterable<IntCoordinates>) this);
		}
		
		@NotNull
		@Override
		public Iterator<IntCoordinates> iterator()
		{
			pointer = new IntCoordinates(lowerBounds);
			pointer.coordinates[pointer.getDimension() - 1]--;
			return this;
		}
		
		@Override
		public boolean hasNext()
		{
			return distanceFromUpper() > pointer.getDimension();
		}
		
		@Override
		public IntCoordinates next()
		{
			pointer = pointer.increaseLastDimension()
			                 .wrap(lowerBounds, upperBounds);
			return pointer;
		}
	}
}
