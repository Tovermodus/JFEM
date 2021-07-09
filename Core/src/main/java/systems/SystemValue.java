package systems;

import basic.PerformanceArguments;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.PressureValue;
import mixed.VelocityValue;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class SystemValue extends DenseVector
{
	int[] starts;
	int[] ends;
	public SystemValue(SystemValue mv)
	{
		super(mv);
		starts = mv.starts.clone();
		ends = mv.ends.clone();
	}
	public SystemValue(int [] ends, double[] values){
		super(values);
		if (PerformanceArguments.getInstance().executeChecks)
			if(ends[ends.length-1] != values.length)
				throw new IllegalArgumentException("last end must be length of values");
		starts = new int[ends.length];
		this.ends = ends.clone();
		for (int i = 0; i < ends.length - 1; i++)
		{
			starts[i+1] = ends[i];
		}
	}
	public SystemValue(int [] ends)
	{
		super(ends[ends.length-1]);
		starts = new int[ends.length];
		this.ends = ends.clone();
		for (int i = 0; i < ends.length - 1; i++)
		{
			starts[i+1] = ends[i];
		}
	}
	public SystemValue()
	{
		this(SystemParameters.getInstance().ends);
	}
	public SystemValue(double [] values)
	{
		this(SystemParameters.getInstance().ends, values);
	}
	public SystemValue(int component, CoordinateVector componentValues)
	{
		this(SystemParameters.getInstance().ends);
		setComponent(componentValues, component);
	}
	public int getNumberOfComponents()
	{
		return ends.length;
	}
	public CoordinateVector getComponent(int component)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (component >= ends.length || component < 0)
				throw new IllegalArgumentException("Component does not exist");
		CoordinateVector ret = new CoordinateVector(ends[component] - starts[component]);
		for (int i = 0; i < ret.getLength(); i++)
		{
			ret.set(entries[starts[component]+i], i);
		}
		return ret;
	}
	public void setComponent(CoordinateVector c, int component)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (component >= ends.length || component < 0)
				throw new IllegalArgumentException("Component does not exist");
			if (ends[component] - starts[component] != c.getLength())
				throw new IllegalArgumentException("Vectors have different size");
		}
		for (int i = 0; i < c.getLength(); i++)
		{
			entries[starts[component]+i] = c.at(i);
		}
	}
	
	@Override
	public String toString()
	{
		return super.toString() +" ends:" + IntStream.of(ends).mapToObj(i-> " "+i).collect(Collectors.joining());
	}
}
