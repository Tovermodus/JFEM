package mixed;

import basic.LagrangeNodeFunctional;
import basic.ScalarFunction;
import basic.VectorFunction;
import com.google.common.collect.Lists;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.SparseMatrix;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.TreeSet;

public class TaylorHoodSpace extends CartesianGridSpace<QkQkFunction>
{
	
	public TaylorHoodSpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
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
		int i = 0;
		for (QkQkFunction shapeFunction : shapeFunctions)
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
				for (TPCell ce : shapeFunction.getCells())
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
			for (int i = 0; i < Math.pow(polynomialDegree + 2, getDimension()) * getDimension(); i++)
			{
				QkQkFunction shapeFunction = new QkQkFunction(new ContinuousTPVectorFunction(cell,
					polynomialDegree + 1, i));
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (TPCell ce : shapeFunction.getCells())
					supportOnCell.put(ce, shapeFunction);
				for (TPFace face : shapeFunction.getFaces())
					supportOnFace.put(face, shapeFunction);
				
			}
		}
	}
	
	public void setVelocityBoundaryValues(VectorFunction boundaryValues)
	{
		setVelocityBoundaryValues(getSystemMatrix());
		setVelocityBoundaryValues(boundaryValues, getRhs());
	}
	public void setVelocityBoundaryValues(SparseMatrix s)
	{
		setVelocityBoundaryValues(ScalarFunction.constantFunction(1),s);
	}
	public void setVelocityBoundaryValues(VectorFunction boundaryValues, DenseVector d)
	{
		setVelocityBoundaryValues(boundaryValues, ScalarFunction.constantFunction(1),d);
	}
	public void setVelocityBoundaryValues(ScalarFunction indicatorFunction, SparseMatrix s)
	{
		List<List<TPFace>> smallerList = Lists.partition(getFaces(), getFaces().size() / 12 + 1);
		smallerList.stream().parallel().forEach(smallList ->
		{
			for (TPFace F : smallList)
			{
				if (F.isBoundaryFace())
				{
					if(TPFaceIntegral.integrateNonTensorProduct(indicatorFunction::value,
						F,QuadratureRule1D.Gauss5) > 0)
					{
						for (MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction,
							ContinuousTPVectorFunction> shapeFunction :
							getShapeFunctionsWithSupportOnFace(F))
						{
							if (shapeFunction.hasVelocityFunction())
							{
								if (F.isOnFace(shapeFunction.getVelocityShapeFunction().getNodeFunctionalPoint()))
								{
									int shapeFunctionIndex = shapeFunction.getGlobalIndex();
									s.deleteLine(shapeFunctionIndex);
									s.set(1, shapeFunctionIndex, shapeFunctionIndex);
								}
							}
						}
					}
				}
			}
		});
	}
	public void setVelocityBoundaryValues(VectorFunction boundaryValues,
	                                      ScalarFunction indicatorFunction, DenseVector d)
	{
		MixedFunction boundaryMixed = new MixedFunction(boundaryValues);
		List<List<TPFace>> smallerList = Lists.partition(getFaces(), getFaces().size() / 12 + 1);
		smallerList.stream().parallel().forEach(smallList ->
		{
			for (TPFace F : smallList)
			{
				if (F.isBoundaryFace())
				{
					if(TPFaceIntegral.integrateNonTensorProduct(indicatorFunction::value,
						F, QuadratureRule1D.Gauss5) > 0)
					{
						for (MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction,
							ContinuousTPVectorFunction> shapeFunction :
							getShapeFunctionsWithSupportOnFace(F))
						{
							if (shapeFunction.hasVelocityFunction())
							{
								double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryMixed);
								if (F.isOnFace(shapeFunction.getVelocityShapeFunction().getNodeFunctionalPoint()))
								{
									int shapeFunctionIndex = shapeFunction.getGlobalIndex();
									d.set(nodeValue, shapeFunctionIndex);
								}
							}
						}
					}
				}
			}
		});
	}
	
	public void setPressureBoundaryValues(ScalarFunction boundaryValues)
	{
		MixedFunction boundaryMixed = new MixedFunction(boundaryValues);
		int progress = 0;
		for (TPFace face : getFaces())
		{
			System.out.println( (100 * progress / getFaces().size()) + "%");
			progress++;
			if (face.isBoundaryFace())
			{
				for (MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction,
					ContinuousTPVectorFunction> shapeFunction :
					getShapeFunctionsWithSupportOnFace(face))
				{
					if (shapeFunction.hasPressureFunction())
					{
						double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryMixed);
						if (nodeValue != 0 || face.isOnFace(((LagrangeNodeFunctional) shapeFunction.getPressureShapeFunction().getNodeFunctional()).getPoint()))
						{
							int shapeFunctionIndex = shapeFunction.getGlobalIndex();
							for (TPCell cell : shapeFunction.getCells())
								for (MixedShapeFunction<TPCell, TPFace,
									ContinuousTPShapeFunction,
									ContinuousTPVectorFunction> sameSupportFunction :
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
	
}
