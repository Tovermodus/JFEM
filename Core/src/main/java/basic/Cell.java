package basic;

import linalg.CoordinateVector;

import java.util.List;
import java.util.Set;

public interface Cell<CT extends Cell<CT,FT>, FT extends Face<CT,FT>> extends Comparable<CT>
{
	int  getDimension();
	
	
	Set<FT> getFaces();
	
	boolean isRefined();
	
	void setRefined(boolean refined);
	
	void addFace(FT face);
	
	boolean isInCell(CoordinateVector pos);
	CoordinateVector center();
	VectorFunction getOuterNormal(FT face);
	
	default List<CT> refine(List<FT> refinedFaces)
	{
		throw new UnsupportedOperationException();
	}

}
