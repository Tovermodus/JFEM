package systems;

import basic.PerformanceArguments;

import java.util.Objects;

public class SystemParameters
{
	private static SystemParameters INSTANCE;
	final int[] ends;
	final int[] starts;
	
	private SystemParameters(int[] ends)
	{
		starts = new int[ends.length];
		this.ends = ends.clone();
		for (int i = 0; i < ends.length - 1; i++)
		{
			starts[i + 1] = ends[i];
		}
	}
	
	public static void createInstance(int[] ends)
	{
		if (INSTANCE == null)
			INSTANCE = new SystemParameters(ends);
		else
			throw new IllegalStateException("Parameters already set");
	}
	
	public static void deleteInstance()
	{
		INSTANCE = null;
	}
	
	public static SystemParameters getInstance()
	{
		if (INSTANCE == null)
			throw new IllegalStateException("No Parameters set");
		return INSTANCE;
	}
	
	
}
