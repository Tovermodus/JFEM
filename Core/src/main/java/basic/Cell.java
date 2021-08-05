package basic;

import com.google.common.collect.ImmutableSet;
import linalg.CoordinateVector;

import java.util.List;

public interface Cell<CT extends Cell<CT,FT>, FT extends Face<CT,FT>> extends Comparable<CT>
{
	int  getDimension();
	
	
	ImmutableSet<FT> getFaces();
	
	
	boolean isInCell(CoordinateVector pos);
	CoordinateVector center();
	VectorFunction getOuterNormal(FT face);
	
	default List<CT> refine(List<FT> refinedFaces)
	{
		throw new UnsupportedOperationException();
	}
	
}
