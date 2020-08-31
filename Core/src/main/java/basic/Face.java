package basic;

import com.google.common.collect.Multimap;
import linalg.CoordinateVector;
import linalg.Vector;

import java.util.List;
import java.util.Set;

public interface Face<CT extends Cell<CT,FT,ST>,FT extends Face<CT,FT,ST>, ST extends ShapeFunction<CT,FT,ST,?,?,?>> extends Comparable<FT>
{
	
	int  getDimension();
	
	Set<CT> getCells();
	
	Set<ST> getShapeFunctions();
	
	void setBoundaryFace(boolean boundaryFace);
	boolean isBoundaryFace();
	
	void addCell(CT cell);
	void addShapeFunction(ST shapeFunction);
	
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
