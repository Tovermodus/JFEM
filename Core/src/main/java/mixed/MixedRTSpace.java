package mixed;

import basic.LagrangeNodeFunctional;
import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;

public class MixedRTSpace extends CartesianGridSpace<RTMixedFunction>
{
	public MixedRTSpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                    List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	@Override
	public void assembleFunctions(int polynomialDegree)
	{
		shapeFunctions = new TreeSet<>();
		assemblePressureFunctions(polynomialDegree);
		assembleVelocityFunctions(polynomialDegree);
	}
	private void assemblePressureFunctions(int polynomialDegree)
	{
		
		for (TPCell cell : getCells())
		{
			for (int i = 0; i < Math.pow(polynomialDegree + 1, getDimension()); i++)
			{
				RTMixedFunction shapeFunction = new RTMixedFunction(new ContinuousTPShapeFunction(cell,
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
			for (int i = 0; i < Math.pow(polynomialDegree + 1, getDimension()-1)*(polynomialDegree+2) * getDimension(); i++)
			{
				RTMixedFunction shapeFunction = new RTMixedFunction(new RTShapeFunction(cell,
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
		int progress = 0;
		for (TPFace face : getFaces())
		{
			System.out.println((100*progress/getFaces().size())+"%");
			progress++;
			if (face.isBoundaryFace())
			{
				for (MixedShapeFunction<TPCell,TPFace, ContinuousTPShapeFunction,
					RTShapeFunction> shapeFunction :
					getShapeFunctionsWithSupportOnFace(face))
				{
					if(shapeFunction.hasVelocityFunction())
					{
						double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryMixed);
						if (nodeValue != 0 || face.isOnFace(shapeFunction.getVelocityShapeFunction().getNodeFunctionalPoint()))
						{
							int shapeFunctionIndex = shapeFunction.getGlobalIndex();
							for (TPCell cell : shapeFunction.getCells())
								for (MixedShapeFunction<TPCell,TPFace,
									ContinuousTPShapeFunction,
									RTShapeFunction> sameSupportFunction :
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
	public void setPressureBoundaryValues(ScalarFunction boundaryValues)
	{
		MixedFunction boundaryMixed = new MixedFunction(boundaryValues);
		for(MixedShapeFunction<TPCell,TPFace, ContinuousTPShapeFunction,RTShapeFunction> shapeFunction :
			getShapeFunctions().values())
		{
			if(!shapeFunction.hasPressureFunction())
				continue;
			boolean boundaryShapeFunction = false;
			for(TPFace f: shapeFunction.getFaces())
			{
				if(f.isBoundaryFace() && f.isOnFace(((LagrangeNodeFunctional)shapeFunction.getPressureShapeFunction().getNodeFunctional()).getPoint()))
				{
					boundaryShapeFunction = true;
					break;
				}
			}
			if(boundaryShapeFunction)
			{
				double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryMixed);
				int shapeFunctionIndex = shapeFunction.getGlobalIndex();
				for (TPCell cell : shapeFunction.getCells())
					for (MixedShapeFunction<TPCell,TPFace, ContinuousTPShapeFunction,
						RTShapeFunction> sameSupportFunction :
						getShapeFunctionsWithSupportOnCell(cell))
						systemMatrix.set(0, shapeFunctionIndex,
							sameSupportFunction.getGlobalIndex());
				getSystemMatrix().set(1, shapeFunctionIndex, shapeFunctionIndex);
				getRhs().set(nodeValue, shapeFunctionIndex);
			}
		}
		
	}
}
