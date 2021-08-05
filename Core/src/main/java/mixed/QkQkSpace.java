package mixed;

import basic.*;
import linalg.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;

public class QkQkSpace extends CartesianGridSpace<QkQkFunction>
{
	
	Set<Integer> boundaryNodes;
	
	public QkQkSpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                 List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
		boundaryNodes = new HashSet<>();
	}
	
	@Override
	public Set<Integer> getFixedNodeIndices()
	{
		return boundaryNodes;
	}
	
	@Override
	public void assembleFunctions(int polynomialDegree)
	{
		shapeFunctions = new TreeSet<>();
		assemblePressureFunctions(polynomialDegree);
		assembleVelocityFunctions(polynomialDegree);
		int i = 0;
		for(QkQkFunction shapeFunction:shapeFunctions)
			shapeFunction.setGlobalIndex(i++);
	}
	private void assemblePressureFunctions(int polynomialDegree)
	{
		
		for (TPCell cell : getCells())
		{
			for (int i = 0; i < Math.pow(polynomialDegree + 1, getDimension()); i++)
			{
				QkQkFunction shapeFunction = new QkQkFunction(new ContinuousTPShapeFunction(cell,
					polynomialDegree, i));
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for(TPCell ce: shapeFunction.getCells())
					supportOnCell.put(ce, shapeFunction);
				for (TPFace face : shapeFunction.getFaces())
					supportOnFace.put(face, shapeFunction);
			}
		}
	}
	private void assembleVelocityFunctions(int polynomialDegree)
	{
		
		for (TPCell cell : getCells())
		{
			for (int i = 0; i < Math.pow(polynomialDegree + 1, getDimension()) * getDimension(); i++)
			{
				QkQkFunction shapeFunction = new QkQkFunction(new ContinuousTPVectorFunction(cell,
					polynomialDegree, i));
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for(TPCell ce: shapeFunction.getCells())
					supportOnCell.put(ce, shapeFunction);
				for (TPFace face : shapeFunction.getFaces())
					supportOnFace.put(face, shapeFunction);
				
			}
		}
	}
	
	public void setVelocityBoundaryValues(VectorFunction boundaryValues)
	{
		MixedFunction boundaryMixed = new MixedFunction(boundaryValues);
		for (TPFace F : getBoundaryFaces())
		{
			for (MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction,
				ContinuousTPVectorFunction> shapeFunction :
				getShapeFunctionsWithSupportOnFace(F))
			{
				if (shapeFunction.hasVelocityFunction())
				{
					double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryMixed);
					if (nodeValue != 0 || F.isOnFace(shapeFunction.getVelocityShapeFunction().getNodeFunctionalPoint()))
					{
						int shapeFunctionIndex = shapeFunction.getGlobalIndex();
						boundaryNodes.add(shapeFunctionIndex);
						getSystemMatrix().deleteLine(shapeFunctionIndex);
						getSystemMatrix().set(1, shapeFunctionIndex, shapeFunctionIndex);
						getRhs().set(nodeValue, shapeFunctionIndex);
					}
				}
			}
			
		}
	}
	
	public void setPressureBoundaryValues(ScalarFunction boundaryValues)
	{
		MixedFunction boundaryMixed = new MixedFunction(boundaryValues);
		for (TPFace F : getBoundaryFaces())
		{
			for (MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction,
				ContinuousTPVectorFunction> shapeFunction :
				getShapeFunctionsWithSupportOnFace(F))
			{
				if (shapeFunction.hasPressureFunction())
				{
					double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryMixed);
					if (nodeValue != 0 || F.isOnFace(((LagrangeNodeFunctional) shapeFunction.getPressureShapeFunction().getNodeFunctional()).getPoint()))
					{
						int shapeFunctionIndex = shapeFunction.getGlobalIndex();
						boundaryNodes.add(shapeFunctionIndex);
						getSystemMatrix().deleteLine(shapeFunctionIndex);
						getSystemMatrix().set(1, shapeFunctionIndex, shapeFunctionIndex);
						getRhs().set(nodeValue, shapeFunctionIndex);
					}
				}
			}
			
		}
	}
}
