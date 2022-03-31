package tensorproduct;

import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.TreeSet;

public class RTFESpace
	extends CartesianGridSpace<RTShapeFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>
{
	public RTFESpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                 final List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	@Override
	public void assembleFunctions(final int polynomialDegree)
	{
		
		shapeFunctions = new TreeSet<>();
		for (final TPCell cell : getCells())
		{
			for (int i = 0; i < Math.pow(polynomialDegree + 1,
			                             getDimension() - 1) * (polynomialDegree + 2) * getDimension(); i++)
			{
				final RTShapeFunction shapeFunction = new RTShapeFunction(cell, polynomialDegree, i);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (final TPCell supportCell : shapeFunction.getCells())
					addFunctionToCell(shapeFunction, supportCell);
				for (final TPFace supportFace : shapeFunction.getFaces())
					getFaceSupportMapping().put(supportFace, shapeFunction);
			}
		}
	}
	
	
	/*public void setBoundaryValues(VectorFunction boundaryValues)
	{
		
		/*int progress = 0;
		for (TPFace face : getFaces())
		{
			System.out.println((100*progress/getFaces().size())+"%");
			progress++;
			if (face.isBoundaryFace())
			{
				for (RTShapeFunction shapeFunction : getShapeFunctionsWithSupportOnFace(face))
				{
					double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryValues);
					if (nodeValue != 0 || face.isOnFace(shapeFunction.getNodeFunctionalPoint()))
					{
						int shapeFunctionIndex = shapeFunction.getGlobalIndex();
						for (TPCell cell : shapeFunction.getCells())
							for (RTShapeFunction sameSupportFunction :
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
