package linalg;

import com.google.common.primitives.Ints;

import java.util.Comparator;
import java.util.List;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.function.ToLongFunction;

public class CoordinateComparator implements Comparator<List<Integer>>
{
	public static int comp(int[] o1, int[] o2)
	{
		if(o1.length != o2.length)
			throw new IllegalArgumentException("coordinates have different length");
		for(int i = 0; i < o1.length; i++)
		{
			if(o1[i] > o2[i])
				return 1;
			if(o1[i] < o2[i])
				return -1;
		}
		return 0;
	}
	public static int comp(double[] o1, double[] o2)
	{
		if(o1.length != o2.length)
			throw new IllegalArgumentException("coordinates have different length");
		for(int i = 0; i < o1.length; i++)
		{
			if(o1[i] > o2[i])
				return 1;
			if(o1[i] < o2[i])
				return -1;
		}
		return 0;
	}
	@Override
	public int compare(List<Integer> o1, List<Integer> o2)
	{
		return comp(Ints.toArray(o1),Ints.toArray(o2));
	}
}
