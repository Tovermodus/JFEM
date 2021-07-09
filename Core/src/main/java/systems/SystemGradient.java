package systems;

import basic.PerformanceArguments;
import linalg.*;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class SystemGradient extends DenseMatrix
{
	int[] starts;
	int[] ends;
	protected SystemGradient(SystemGradient mv)
	{
		super(mv);
		starts = mv.starts.clone();
		ends = mv.ends.clone();
	}
	protected SystemGradient(int [] ends)
	{
		super(ends[ends.length-1], ends[ends.length-1]);
		starts = new int[ends.length];
		this.ends = ends.clone();
		for (int i = 0; i < ends.length - 1; i++)
		{
			starts[i+1] = ends[i];
		}
	}
	protected SystemGradient()
	{
		this(SystemParameters.getInstance().ends);
	}
	
	public int getNumberOfComponents()
	{
		return ends.length;
	}
	public CoordinateMatrix getComponent(int componentY, int componentX)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (componentX >= ends.length || componentX < 0)
				throw new IllegalArgumentException("Component does not exist");
			if (componentY >= ends.length || componentY < 0)
				throw new IllegalArgumentException("Component does not exist");
		}
		CoordinateMatrix ret = new CoordinateMatrix(ends[componentY] - starts[componentY],
			ends[componentX] - starts[componentX]);
		for (int i = 0; i < ret.getRows(); i++)
		{
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(entries[starts[componentY]+i][starts[componentX]+j], i, j);
		}
		return ret;
	}
	public void setComponent(CoordinateMatrix c, int componentY, int componentX)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (componentX >= ends.length || componentX < 0)
				throw new IllegalArgumentException("Component does not exist");
			if (componentY >= ends.length || componentY < 0)
				throw new IllegalArgumentException("Component does not exist");
			if (!new IntCoordinates(ends[componentY] - starts[componentY],
				ends[componentX] - starts[componentX]).equals(c.getShape()))
				throw new IllegalArgumentException("Vectors have different size");
		}
		for (int i = 0; i < c.getRows(); i++)
		{
			for (int j = 0; j < c.getCols(); j++)
				entries[starts[componentY]+i][starts[componentX]+j] = c.at(i, j);
		}
	}
	public void setComponent(double c, int componentY, int componentX)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (componentX >= ends.length || componentX < 0)
				throw new IllegalArgumentException("Component does not exist");
			if (componentY >= ends.length || componentY < 0)
				throw new IllegalArgumentException("Component does not exist");
			if (!new IntCoordinates(ends[componentY] - starts[componentY],
				ends[componentX] - starts[componentX]).equals(new IntCoordinates(1,1)))
				throw new IllegalArgumentException("Vectors have different size");
		}
		entries[starts[componentY]][starts[componentX]] = c;
	}
	
	
	
	@Override
	public String toString()
	{
		return super.toString() +" ends:" + IntStream.of(ends).mapToObj(i-> " "+i).collect(Collectors.joining());
	}
}
