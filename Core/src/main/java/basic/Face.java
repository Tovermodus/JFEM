package basic;

import com.google.common.collect.Multimap;
import linalg.CoordinateVector;
import linalg.Vector;

import java.util.List;
import java.util.Set;

public interface Face<CT extends Cell<CT,FT, ET>,FT extends Face<CT,FT, ET>, ET extends Edge<CT,FT,ET>> extends Comparable<FT>
{
	
	int  getDimension();
	
	Set<CT> getCells();
	
	
	void setBoundaryFace(boolean boundaryFace);
	boolean isBoundaryFace();
	
	void addCell(CT cell);
	
	VectorFunction getNormal();
	
	CT getNormalDownstreamCell(CoordinateVector pos);
	CT getNormalUpstreamCell(CoordinateVector pos);
	
	CoordinateVector center();
	
	boolean isOnFace(CoordinateVector pos);
	
	default List<FT> refine(Multimap<CT, CT> cellMap)
	{
		throw new UnsupportedOperationException();
	}
	
}
