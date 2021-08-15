package distorted;

import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.TPShapeFunction;

import java.util.HashSet;
import java.util.TreeSet;

/*
Workflow of this finite element is: first generate 1D cells, then
 */
public class DistortedSpace extends CircleGridSpace<DistortedShapeFunction, Double, CoordinateVector, CoordinateMatrix>
{
	public DistortedSpace(final CoordinateVector center, final double radius, final int refinements)
	{
		super(center, radius, refinements);
	}
	
	@Override
	public void assembleFunctions(final int polynomialDegree)
	{
		shapeFunctions = new TreeSet<>();
		for (final DistortedCell cell : grid.cells)
		{
			for (int localIndex = 0; localIndex < TPShapeFunction.functionsPerCell(polynomialDegree,
			                                                                       getDimension()); localIndex++)
			{
				final DistortedShapeFunction function = new DistortedShapeFunction(cell,
				                                                                   polynomialDegree,
				                                                                   localIndex);
				function.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(function);
				for (final DistortedCell c : function.getCells())
					supportOnCell.put(c, function);
				for (final DistortedFace f : function.getFaces())
					supportOnFace.put(f, function);
			}
		}
	}
}
