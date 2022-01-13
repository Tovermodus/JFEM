package tensorproduct;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class ContinuousTPFESpace
	extends CartesianGridSpace<ContinuousTPShapeFunction, Double, CoordinateVector, CoordinateMatrix>
{
	
	public ContinuousTPFESpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                           final List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	public ContinuousTPFESpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
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
				final ContinuousTPShapeFunction shapeFunction = new ContinuousTPShapeFunction(cell,
				                                                                              polynomialDegree,
				                                                                              i);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (final TPCell supportCell : shapeFunction.getCells())
					supportOnCell.put(supportCell, shapeFunction);
				for (final TPFace supportFace : shapeFunction.getFaces())
					supportOnFace.put(supportFace, shapeFunction);
			}
		if (shapeFunctions.size() != Arrays.stream(grid.cellsPerDimension.asArray())
		                                   .map(cellsPerDimension -> cellsPerDimension * polynomialDegree + 1)
		                                   .reduce(1, Math::multiplyExact))
			throw new IllegalStateException("Identification did not work");
	}

//	public void setBoundaryValues(ScalarFunction boundaryValues)
//	{
//		for (TPFace face : getFaces())
//		{
//			if (face.isBoundaryFace())
//			{
//				for (ContinuousTPShapeFunction shapeFunction : getShapeFunctionsWithSupportOnFace(face))
//				{
//					double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryValues);
//					if (nodeValue != 0 || face.isOnFace(((LagrangeNodeFunctional)shapeFunction.getNodeFunctional()).getPoint()))
//					{
//						int shapeFunctionIndex = shapeFunction.getGlobalIndex();
//						for (TPCell cell : shapeFunction.getCells())
//							for (ContinuousTPShapeFunction sameSupportFunction :
//								getShapeFunctionsWithSupportOnCell(cell))
//								systemMatrix.set(0, shapeFunctionIndex,
//									sameSupportFunction.getGlobalIndex());
//						getSystemMatrix().set(1, shapeFunctionIndex, shapeFunctionIndex);
//						getRhs().set(nodeValue, shapeFunctionIndex);
//					}
//				}
//			}
//		}
//	}
}
