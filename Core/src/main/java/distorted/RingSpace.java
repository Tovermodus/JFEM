package distorted;

import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.TPShapeFunction;

import java.util.TreeSet;

/*
Workflow of this finite element is: first generate 1D cells, then
 */
public class RingSpace
	extends RingGridSpace<DistortedShapeFunction, Double, CoordinateVector, CoordinateMatrix>
{
	public RingSpace(final CoordinateVector center,
	                 final double innerRadius,
	                 final double outerRadius,
	                 final int refinements)
	{
		super(center, innerRadius, outerRadius, refinements);
	}
	
	@Override
	public void assembleFunctions(final int polynomialDegree)
	{
		shapeFunctions = new TreeSet<>();
		for (final DistortedCell cell : grid.getCells())
		{
			for (int localIndex = 0; localIndex < TPShapeFunction.functionsPerCell(polynomialDegree,
			                                                                       getDimension()); localIndex++)
			{
				final DistortedShapeFunction function = new DistortedShapeFunction(cell,
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
