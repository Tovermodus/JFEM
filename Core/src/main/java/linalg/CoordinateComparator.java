package linalg;

import com.google.common.primitives.Ints;

import java.util.Comparator;
import java.util.List;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.function.ToLongFunction;

public class CoordinateComparator
{
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
	public static int comp(CoordinateVector o1, CoordinateVector o2)
	{
		if(o1.getLength() != o2.getLength())
			throw new IllegalArgumentException("coordinates have different length");
		for(int i = 0; i < o1.getLength(); i++)
		{
			if(o1.at(i) > o2.at(i)+1e-15)
				return 1;
			if(o1.at(i) < o2.at(i)-1e-15)
				return -1;
		}
		return 0;
	}
}
