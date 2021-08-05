package tensorproduct;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.HashSet;
import java.util.List;

/*
Workflow of this finite element is: first generate 1D cells, then
 */
public class TPFESpace extends CartesianGridSpace<TPShapeFunction>
{
	public TPFESpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                 List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	@Override
	public void assembleFunctions(int polynomialDegree)
	{
		shapeFunctions = new HashSet<>();
		for(TPCell cell: grid.cells)
		{
			for(int localIndex = 0; localIndex < TPShapeFunction.functionsPerCell(polynomialDegree, getDimension()); localIndex++)
			{
				TPShapeFunction function = new TPShapeFunction(cell, polynomialDegree, localIndex);
				function.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(function);
				supportOnCell.put(cell, function);
				for(TPFace f:cell.getFaces())
					supportOnFace.put(f, function);
			}
		}
	}

//	/*
//	TODO: Provide interface for fast TPIntegrals
//	 */
//
}
