package tensorproduct;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class CTPFESpace
	extends CartesianGridSpace<CTPShapeFunction, Double, CoordinateVector, CoordinateMatrix>
{
	
	public CTPFESpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                  final List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	public CTPFESpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                  final IntCoordinates cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	@Override
	public void assembleFunctions(final int polynomialDegree)
	{
		shapeFunctions = new HashSet<>();
		for (final TPCell cell : grid.cells)
			for (int i = 0; i < Math.pow(polynomialDegree + 1, getDimension()); i++)
			{
				final CTPShapeFunction shapeFunction = new CTPShapeFunction(cell,
				                                                            polynomialDegree,
				                                                            i);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (final TPCell supportCell : shapeFunction.getCells())
					addFunctionToCell(shapeFunction, supportCell);
				for (final TPFace supportFace : shapeFunction.getFaces())
					getFaceSupportMapping().put(supportFace, shapeFunction);
			}
		if (shapeFunctions.size() != Arrays.stream(grid.cellsPerDimension.asArray())
		                                   .map(cellsPerDimension -> cellsPerDimension * polynomialDegree + 1)
		                                   .reduce(1, Math::multiplyExact))
			throw new IllegalStateException("Identification did not work");
	}
}
