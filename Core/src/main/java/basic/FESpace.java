package basic;

import com.google.common.collect.Multimap;
import linalg.CoordinateVector;

import java.util.*;

public interface FESpace<CT extends Cell<CT,FT>, FT extends  Face<CT,FT>,
	ST extends ShapeFunction<CT,FT,?,?,?>, FST extends FESpace<CT,FT, ST,FST>>
{
	int getDimension();
	List<CT> getCells();
	Map<Integer, ST> getShapeFunctions();
	List<FT> getFaces();
	default List<FT> getBoundaryFaces()
	{
		List<FT> boundaryFaces = new ArrayList<>();
		for (FT F : getFaces())
			if (F.isBoundaryFace())
				boundaryFaces.add(F);
		return boundaryFaces;
	}
	Collection<ST> getShapeFunctionsWithSupportOnCell(CT cell);
	Collection<ST> getShapeFunctionsWithSupportOnFace(FT face);
	
	List<CoordinateVector> generatePlotPoints(int resolution);
	default FST refine(Multimap<CT, CT> cellRefinedCellMapping,
	           Multimap<FT, FT> faceRefinedFaceMapping)
	{
		throw new UnsupportedOperationException();
	}
}
