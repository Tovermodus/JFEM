package distorted;

import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;
import tensorproduct.TPShapeFunction;

import java.util.TreeSet;

/*
Workflow of this finite element is: first generate 1D cells, then
 */
public class CircleVectorSpace
	extends CircleGridSpace<DistortedVectorShapeFunction, CoordinateVector,
	CoordinateMatrix, CoordinateTensor>
{
	public CircleVectorSpace(final CoordinateVector center, final double radius, final int refinements)
	{
		super(center, radius, refinements);
	}
	
	@Override
	public void assembleFunctions(final int polynomialDegree)
	{
		shapeFunctions = new TreeSet<>();
		for (final DistortedCell cell : grid.cells)
		{
			for (int localIndex = 0; localIndex < getDimension() * TPShapeFunction.functionsPerCell(
				polynomialDegree,
				getDimension()); localIndex++)
			{
				final DistortedVectorShapeFunction function = new DistortedVectorShapeFunction(cell,
				                                                                               polynomialDegree,
				                                                                               localIndex);
				if (shapeFunctions.add(function))
					function.setGlobalIndex(shapeFunctions.size() - 1);
				
				for (final DistortedCell c : function.getCells())
				{
					addFunctionToCell(function, c);
				}
				for (final DistortedFace f : function.getFaces())
					addFunctionToFace(function, f);
			}
		}
	}
}
