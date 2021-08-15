package basic;

import com.google.common.collect.ImmutableSet;
import linalg.CoordinateVector;

import java.util.List;

public interface Cell<CT extends Cell<CT, FT>, FT extends Face<CT, FT>> extends Comparable<CT>
{
	int getDimension();
	
	ImmutableSet<FT> getFaces();
	
	boolean isInCell(CoordinateVector pos);
	
	CoordinateVector center();
	
	VectorFunction getOuterNormal(FT face);
	
	double diam();
	
	default ScalarFunction indicatorFunction()
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return getDimension();
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				if (isInCell(pos))
					return 1.0;
				return 0.0;
			}
		};
	}
	
	default List<CT> refine(final List<FT> refinedFaces)
	{
		throw new UnsupportedOperationException();
	}
}
