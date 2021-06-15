package linalg;

import org.jetbrains.annotations.NotNull;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

public class IntCoordinates implements Comparable<IntCoordinates>, Cloneable
{
	final private int [] coordinates;
	public IntCoordinates(int ... coordinates) {
		this.coordinates = coordinates;
	}
	public IntCoordinates(IntCoordinates source) {
		this.coordinates = source.coordinates.clone();
	}
	public int [] getCoordinates() {
		return coordinates;
	}
	public int getDimension() {
		return coordinates.length;
	}
	private IntCoordinates increaseLastDimension() {
		int [] ret = coordinates;
		ret[getDimension() - 1]++;
		return new IntCoordinates(ret);
	}
	private IntCoordinates wrap(IntCoordinates lowerBounds, IntCoordinates upperBounds) {
		for(int d = getDimension() - 1; d > 0; d--)
			if(coordinates[d] == upperBounds.coordinates[d])
			{
				coordinates[d] = lowerBounds.coordinates[d];
				coordinates[d-1]++;
			}
		return this;
	}
	@Override
	public int hashCode()
	{
		int ret = 0;
		for(int c:coordinates)
		{
			ret *= 104729;
			ret += c;
			
		}
		return ret;
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(obj instanceof IntCoordinates)
			return this.compareTo((IntCoordinates) obj) == 0;
		return false;
	}
	
	@Override
	public String toString()
	{
		StringBuilder ret = new StringBuilder();
		for(int c: coordinates)
			ret.append(c).append(" ");
		return ret.toString();
	}
	
	@Override
	public int compareTo(@NotNull IntCoordinates o)
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
	
	public static class Range implements Iterable<IntCoordinates>, Iterator<IntCoordinates>
	{
		final IntCoordinates lowerBounds;
		final IntCoordinates upperBounds;
		IntCoordinates pointer;
		public Range(int[] lowerBounds, int[] upperBounds)
		{
			this(new IntCoordinates(lowerBounds), new IntCoordinates(upperBounds));
		}
		public Range(IntCoordinates lowerBounds, IntCoordinates upperBounds)
		{
			if(lowerBounds.getDimension() != upperBounds.getDimension())
				throw new IllegalArgumentException("Range: bounds must have same length!");
			for (int i = 0; i < lowerBounds.getDimension(); i++)
			{
				if(lowerBounds.getCoordinates()[i] >= upperBounds.getCoordinates()[i])
					throw new IllegalArgumentException("Range: lower index higher than upper!");
			}
			this.lowerBounds = lowerBounds;
			this.upperBounds = upperBounds;
			pointer = new IntCoordinates(lowerBounds);
		}
		
		@NotNull
		@Override
		public Iterator<IntCoordinates> iterator()
		{
			pointer = lowerBounds;
			return this;
		}
		
		@Override
		public boolean hasNext()
		{
			return pointer.increaseLastDimension().wrap(lowerBounds, upperBounds).equals(upperBounds);
		}
		@Override
		public IntCoordinates next()
		{
			return pointer.increaseLastDimension().wrap(lowerBounds, upperBounds);
		}
	}
	
	
}