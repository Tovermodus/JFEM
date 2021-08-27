package systems;

import basic.PerformanceArguments;
import linalg.*;

import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SystemGradient extends DenseMatrix
{
	final int[] starts;
	final int[] ends;
	
	public SystemGradient(final SystemGradient mv)
	{
		super(mv);
		starts = mv.starts;
		ends = mv.ends;
	}
	
	public SystemGradient(final int[] ends, final int d)
	{
		super(ends[ends.length - 1], d);
		starts = new int[ends.length];
		this.ends = ends.clone();
		for (int i = 0; i < ends.length - 1; i++)
		{
			starts[i + 1] = ends[i];
		}
	}
	
	public SystemGradient(final int d)
	{
		super(SystemParameters.getInstance().ends[SystemParameters.getInstance().ends.length - 1], d);
		starts = SystemParameters.getInstance().starts;
		ends = SystemParameters.getInstance().ends;
	}
	
	public SystemGradient(final DenseMatrix d, final boolean wrap)
	{
		super(d, wrap);
		starts = SystemParameters.getInstance().starts;
		ends = SystemParameters.getInstance().ends;
	}
	
	public int getNumberOfComponents()
	{
		return ends.length;
	}
	
	public CoordinateMatrix getComponent(final int component)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (component >= ends.length || component < 0)
				throw new IllegalArgumentException("Component does not exist");
		}
		final CoordinateDenseMatrix ret = new CoordinateDenseMatrix(ends[component] - starts[component], getCols());
		for (int i = 0; i < ret.getRows(); i++)
		{
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(entries[starts[component] + i][j], i, j);
		}
		return ret;
	}
	
	public void setComponent(final CoordinateMatrix c, final int component)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (component >= ends.length || component < 0)
				throw new IllegalArgumentException("Component does not exist");
			if (!new IntCoordinates(ends[component] - starts[component], getCols()).equals(c.getShape()))
				throw new IllegalArgumentException("Vectors have different size");
		}
		for (int i = 0; i < c.getRows(); i++)
		{
			for (int j = 0; j < c.getCols(); j++)
				entries[starts[component] + i][j] = c.at(i, j);
		}
	}
	
	public void setComponent(final CoordinateVector c, final int component)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (component >= ends.length || component < 0)
				throw new IllegalArgumentException("Component does not exist");
			if (ends[component] - starts[component] != 1)
				throw new IllegalArgumentException("Vectors have different size");
		}
		for (int i = 0; i < c.getLength(); i++)
		{
			entries[starts[component]][i] = c.at(i);
		}
	}
	
	@Override
	public String toString()
	{
		return super.toString() + " ends:" + IntStream
			.of(ends)
			.mapToObj(i -> " " + i)
			.collect(Collectors.joining());
	}
}
