package basic;

import com.google.common.collect.Multimap;
import linalg.CoordinateVector;
import tensorproduct.TPCell;
import tensorproduct.TPFace;

import java.util.List;
import java.util.Set;

public interface Edge<CT extends Cell<CT,FT,ET>,FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>> extends Comparable<ET>
{
	
	void addCell(TPCell cell);
	
	void addFace(TPFace cell);
	
	int  getDimension();
	
	Set<CT> getCells();
	Set<FT> getFaces();
	
	void setBoundaryFace(boolean boundaryFace);
	boolean isBoundaryEdge();
	
	VectorFunction getTangent();
	
	CoordinateVector center();
	
	boolean isOnEdge(CoordinateVector pos);
	
	default List<ET> refine(Multimap<CT, CT> cellMap)
	{
		throw new UnsupportedOperationException();
	}
	
}
