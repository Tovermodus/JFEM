package basic;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

import java.util.List;

public interface Face<CT extends Cell<CT,FT>,FT extends Face<CT,FT>> extends Comparable<FT>
{
	
	int  getDimension();
	
	ImmutableSet<CT> getCells();
	
	boolean isBoundaryFace();
	
	VectorFunction getNormal();
	
	CT getNormalDownstreamCell();
	CT getNormalUpstreamCell();
	
	CoordinateVector center();
	
	boolean isOnFace(CoordinateVector pos);
	
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
			public Double value(CoordinateVector pos)
			{
				if(isOnFace(pos))
					return 1.0;
				return 0.0;
			}
		};
	}
	
	default List<FT> refine(Multimap<CT, CT> cellMap)
	{
		throw new UnsupportedOperationException();
	}
	
}
