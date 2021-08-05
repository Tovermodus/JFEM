package tensorproduct;

import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

public class ContinuousTPFEVectorSpace extends CartesianGridSpace<ContinuousTPVectorFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>
{
	
	public ContinuousTPFEVectorSpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                                 List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	@Override
	public void assembleFunctions(int polynomialDegree)
	{
		
		shapeFunctions = new TreeSet<>();
		for (TPCell cell : grid.cells)
		{
			for (int i = 0; i < Math.pow(polynomialDegree + 1, getDimension()) * getDimension(); i++)
			{
				ContinuousTPVectorFunction shapeFunction = new ContinuousTPVectorFunction(cell, polynomialDegree, i);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (TPCell supportCell : shapeFunction.getCells())
					supportOnCell.put(supportCell, shapeFunction);
				for (TPFace supportFace : shapeFunction.getFaces())
					supportOnFace.put(supportFace, shapeFunction);
			}
		}
		if (shapeFunctions.size() != getDimension()*Arrays.stream(grid.cellsPerDimension.asArray())
			.map(cellsPerDimension -> cellsPerDimension * polynomialDegree + 1)
			.reduce(1, Math::multiplyExact))
			throw new IllegalStateException("Identification did not work");
	}
	
	/*
	public void setBoundaryValues(VectorFunction boundaryValues)
	{
		int progress = 0;
		for (TPFace face : getFaces())
		{
			System.out.println((100*progress/getFaces().size())+"%");
			progress++;
			if (face.isBoundaryFace())
			{
				for (ContinuousTPVectorFunction shapeFunction : getShapeFunctionsWithSupportOnFace(face))
				{
					double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryValues);
					if (nodeValue != 0 || face.isOnFace(shapeFunction.getNodeFunctionalPoint()))
					{
						int shapeFunctionIndex = shapeFunction.getGlobalIndex();
						for (TPCell cell : shapeFunction.getCells())
							for (ContinuousTPVectorFunction sameSupportFunction :
								getShapeFunctionsWithSupportOnCell(cell))
								systemMatrix.set(0, shapeFunctionIndex,
									sameSupportFunction.getGlobalIndex());
						getSystemMatrix().set(1, shapeFunctionIndex, shapeFunctionIndex);
						getRhs().set(nodeValue, shapeFunctionIndex);
					}
				}
			}
		}
	}*/
}

