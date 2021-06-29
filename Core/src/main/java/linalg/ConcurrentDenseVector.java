package linalg;


import basic.PerformanceArguments;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

public class ConcurrentDenseVector implements MutableVector
{
	double[] entries;
	BlockingQueue<IntCoordinatesValue> writtenElements;
	ArrayList<BlockingQueue<IntCoordinatesValue>> cacheWrites;
	
	public ConcurrentDenseVector(int size)
	{
		this(new double[size]);
	}
	
	public ConcurrentDenseVector(Vector vect)
	{
		this(vect.getLength());
		for (int i = 0; i < vect.getLength(); i++)
		{
			entries[i] = vect.at(i);
		}
	}
	
	public ConcurrentDenseVector(double[] vect)
	{
		entries = vect.clone();
		writtenElements =
			new ArrayBlockingQueue<>(entries.length/5+100);
		cacheWrites = new ArrayList<>(entries.length/PerformanceArguments.getInstance().cacheSize+1);
		for (int i = 0; i < entries.length/PerformanceArguments.getInstance().cacheSize+1; i++)
		{
			cacheWrites.add(new ArrayBlockingQueue<>(PerformanceArguments.getInstance().cacheSize));
		}
	}
	
	public static ConcurrentDenseVector vectorFromValues(double... values)
	{
		return new ConcurrentDenseVector(values);
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		set(value, new IntCoordinates(coordinates));
	}
	@Override
	public void set(double value, IntCoordinates coordinates)
	{
		try
		{
			writtenElements.put(new IntCoordinatesValue(coordinates, value));
		} catch (InterruptedException e)
		{
			throw new IllegalStateException(" Thread was interrupted while writing to vector"+e);
		}
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public void addInPlace(Tensor other)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public void mulInPlace(double scalar)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public double at(int... coordinates)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int getSparseEntryCount()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public boolean isSparse()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public IntCoordinates getShape()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Vector add(Tensor other)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Vector mul(double scalar)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Matrix outer(Vector other)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public double inner(Vector other)
	{
		return MutableVector.super.inner(other);
	}
	
	@Override
	public void setAll(Vector values, int startingPoint)
	{
		MutableVector.super.setAll(values, startingPoint);
	}
	
	@Override
	public void addAll(Vector values, int startingPoint)
	{
		MutableVector.super.addAll(values, startingPoint);
	}
}
