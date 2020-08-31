package basic;

import com.google.common.collect.Multimap;
import linalg.CoordinateVector;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

public interface FESpace<CT extends Cell<CT,FT>, FT extends  Face<CT,FT>, ST extends ShapeFunction<CT,FT,
	ST,valueT,gradientT,hessianT>,valueT,gradientT,hessianT, FST extends FESpace<CT,FT,ST,valueT,gradientT,
	hessianT,FST>>
{
	List<CT> getCells();
	Map<Integer, ST> getShapeFunctions();
	List<FT> getFaces();
	Set<ST> getShapeFunctionsWithSupportOnCell(CT cell);
	Set<ST> getShapeFunctionsWithSupportOnFace(FT face);
	
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
