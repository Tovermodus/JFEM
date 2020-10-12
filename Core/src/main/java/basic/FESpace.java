package basic;

import com.google.common.collect.Multimap;
import linalg.CoordinateVector;

import java.util.*;

public interface FESpace<CT extends Cell<CT,FT,ET>, FT extends  Face<CT,FT,ET>,
	ET extends Edge<CT,FT,ET>,ST extends ShapeFunction<CT,FT,
	ET,ST,valueT,gradientT,hessianT>,valueT,gradientT,hessianT, FST extends FESpace<CT,FT,ET,ST,valueT,gradientT,
	hessianT,FST>>
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
//	{
//		FESpace ret = new FESpace();
//		ArrayList<Cell> refinedCells = ret.cells;
//		ArrayList<Face> refinedFaces = ret.faces;
//		ArrayList<ScalarShapeFunction> refinedShapeFunctions = ret.shapeFunctions;
//		for(Cell cell:cells)
//		{
//			ArrayList<Cell> refCells = cell.refine(refinedFaces);
//			for(Cell refCell: refCells)
//			{
//				refinedCells.add(refCell);
//				cellRefinedCellMapping.put(cell, refCell);
//			}
//		}
//		for(Face f:faces)
//		{
//			ArrayList<Face> refFaces = f.refine(cellRefinedCellMapping);
//			for(Face refFace:refFaces)
//			{
//				refinedFaces.add(refFace);
//				faceRefinedFaceMapping.put(f,refFace);
//			}
//		}
//		for(Cell refinedCell: refinedCells)
//		{
//			refinedCell.distributeFunctions(refinedShapeFunctions);
//		}
//		return ret;
//
//	}


}
