package basic;

import linalg.CoordinateVector;

import java.util.List;
import java.util.Set;

public interface Cell<CT extends Cell<CT,FT,ST>, FT extends Face<CT,FT,ST>, ST extends ShapeFunction<CT,FT,ST,?,?,?>> extends Comparable<CT>
{
	int  getDimension();
	
	Set<ST> getShapeFunctions();
	
	Set<FT> getFaces();
	
	boolean isRefined();
	
	void setRefined(boolean refined);
	
	void addFace(FT face);
	void addShapeFunction(ST shapeFunction);
	
	boolean isInCell(CoordinateVector pos);
	CoordinateVector center();
	
	default List<CT> refine(List<FT> refinedFaces)
	{
		throw new UnsupportedOperationException();
	}

}
