package tensorproduct;

import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;

public class TPVectorFESpace extends CartesianGridSpace<TPVectorFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>
{
	public TPVectorFESpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                       final List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	@Override
	public void assembleFunctions(final int polynomialDegree)
	{
		shapeFunctions = new HashSet<>();
		for (final TPCell cell : getCells())
		{
			for (int i = 0; i < Math.pow(polynomialDegree + 1, getDimension()) * getDimension(); i++)
			{
				final TPVectorFunction shapeFunction = new TPVectorFunction(cell, polynomialDegree, i,
				                                                            TPShapeFunction.class);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				supportOnCell.put(cell, shapeFunction);
				for (final TPFace face : cell.getFaces())
					supportOnFace.put(face, shapeFunction);
			}
		}
	}
}
