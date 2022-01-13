package tensorproduct;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.TreeSet;

public class ScalarRTFESpace
	extends CartesianGridSpace<RTComponentFunction, Double, CoordinateVector, CoordinateMatrix>
{
	public ScalarRTFESpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                       final List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	@Override
	public void assembleFunctions(final int polynomialDegree)
	{
		shapeFunctions = new TreeSet<>();
		for (final TPCell cell : getCells())
			for (int i = 0; i < Math.pow(polynomialDegree + 1,
			                             getDimension() - 1) * (polynomialDegree + 2); i++)
			{
				final RTComponentFunction shapeFunction = new RTComponentFunction(cell,
				                                                                  polynomialDegree, i, 0);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (final TPCell supportCell : shapeFunction.getCells())
					getCellSupportMapping().put(supportCell, shapeFunction);
				for (final TPFace supportFace : shapeFunction.getFaces())
					getFaceSupportMapping().put(supportFace, shapeFunction);
			}
		if (getDimension() == 2)
			if (shapeFunctions.size() != ((polynomialDegree + 1) * grid.cellsPerDimension.get(0) + 1) * (polynomialDegree + 1) * grid.cellsPerDimension.get(
				1))
				throw new IllegalStateException("Identification did not work ");
		if (getDimension() == 3)
			if (shapeFunctions.size() != ((polynomialDegree + 1) * grid.cellsPerDimension.get(0) + 1)
				* (polynomialDegree + 1) * grid.cellsPerDimension.get(1)
				* (polynomialDegree + 1) * grid.cellsPerDimension.get(2))
				throw new IllegalStateException("Identification did not work");
	}
}
