package linalg;

import org.jetbrains.annotations.NotNull;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

public class IntCoordinates implements Comparable<IntCoordinates>
{
	final private int [] coordinates;
	public IntCoordinates(int ... coordinates) {
		this.coordinates = coordinates;
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
	private IntCoordinates wrap(int[] lowerBounds, int[] upperBounds) {
		for(int d = getDimension() - 1; d > 0; d--)
			if(coordinates[d] == upperBounds[d])
			{
				coordinates[d] = lowerBounds[d];
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
	public class Range implements Iterable<IntCoordinates>, Iterator<IntCoordinates>
	{
		final int[] lowerBounds;
		final int[] upperBounds;
		IntCoordinates pointer;
		public Range(int[] lowerBounds, int[] upperBounds)
		{
			this.lowerBounds = lowerBounds;
			this.upperBounds = upperBounds;
			pointer = new IntCoordinates(lowerBounds);
		}
		
		@NotNull
		@Override
		public Iterator<IntCoordinates> iterator()
		{
			pointer = new IntCoordinates(lowerBounds);
			return this;
		}
		
		@Override
		public boolean hasNext()
		{
			return Arrays.equals(pointer.increaseLastDimension().wrap(lowerBounds, upperBounds).coordinates,
				upperBounds);
		}
		
		@Override
		public IntCoordinates next()
		{
			return pointer.increaseLastDimension().wrap(lowerBounds, upperBounds);
		}
	}
	
	
}