package tensorproduct;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.HashSet;
import java.util.List;

/*
Workflow of this finite element is: first generate 1D cells, then
 */
public class TPFESpace
	extends CartesianGridSpace<TPShapeFunction, Double, CoordinateVector, CoordinateMatrix>
{
	public TPFESpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                 final List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	public TPFESpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                 final IntCoordinates cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	@Override
	public void assembleFunctions(final int polynomialDegree)
	{
		shapeFunctions = new HashSet<>();
		for (final TPCell cell : grid.cells)
		{
			for (int localIndex = 0; localIndex < TPShapeFunction.functionsPerCell(polynomialDegree,
			                                                                       getDimension()); localIndex++)
			{
				final TPShapeFunction function = new TPShapeFunction(cell,
				                                                     polynomialDegree,
				                                                     localIndex);
				function.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(function);
				getCellSupportMapping().put(cell, function);
				for (final TPFace f : cell.getFaces())
					getFaceSupportMapping().put(f, function);
			}
		}
	}

//	/*
//	TODO: Provide interface for fast TPIntegrals
//	 */
//
}
