package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.ScalarFunction;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

public class ContinuousTPFESpace extends CartesianGridSpace<ContinuousTPShapeFunction>
{
	
	public ContinuousTPFESpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                 List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	@Override
	public void assembleFunctions(int polynomialDegree)
	{
		shapeFunctions = new HashSet<>();
		for (TPCell cell : grid.cells)
			for (int i = 0; i < Math.pow(polynomialDegree + 1, getDimension()); i++)
			{
				ContinuousTPShapeFunction shapeFunction = new ContinuousTPShapeFunction(cell,
					polynomialDegree, i);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (TPCell supportCell : shapeFunction.getCells())
					supportOnCell.put(supportCell, shapeFunction);
				for (TPFace supportFace : shapeFunction.getFaces())
					supportOnFace.put(supportFace, shapeFunction);
			}
		if (shapeFunctions.size() != Arrays.stream(grid.cellsPerDimension.asArray())
			.map(cellsPerDimension -> cellsPerDimension * polynomialDegree + 1)
			.reduce(1, Math::multiplyExact))
			throw new IllegalStateException("Identification did not work");
	}
	
	public void setBoundaryValues(ScalarFunction boundaryValues)
	{
		for (TPFace face : getFaces())
		{
			if (face.isBoundaryFace())
			{
				for (ContinuousTPShapeFunction shapeFunction : getShapeFunctionsWithSupportOnFace(face))
				{
					double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryValues);
					if (nodeValue != 0 || face.isOnFace(((LagrangeNodeFunctional)shapeFunction.getNodeFunctional()).getPoint()))
					{
						int shapeFunctionIndex = shapeFunction.getGlobalIndex();
						for (TPCell cell : shapeFunction.getCells())
							for (ContinuousTPShapeFunction sameSupportFunction :
								getShapeFunctionsWithSupportOnCell(cell))
								systemMatrix.set(0, shapeFunctionIndex,
									sameSupportFunction.getGlobalIndex());
						getSystemMatrix().set(1, shapeFunctionIndex, shapeFunctionIndex);
						getRhs().set(nodeValue, shapeFunctionIndex);
					}
				}
			}
		}
	}
	
}

