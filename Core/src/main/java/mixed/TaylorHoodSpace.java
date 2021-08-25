package mixed;

import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.SparseMatrix;
import tensorproduct.CartesianGridSpace;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.TreeSet;

public class TaylorHoodSpace extends CartesianGridSpace<QkQkFunction, MixedValue, MixedGradient, MixedHessian>
{
	
	public TaylorHoodSpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                       final List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	public TaylorHoodSpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                       final IntCoordinates cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	@Override
	public void assembleFunctions(final int polynomialDegree)
	{
		shapeFunctions = new TreeSet<>();
		assemblePressureFunctions(polynomialDegree);
		assembleVelocityFunctions(polynomialDegree);
		int i = 0;
		for (final QkQkFunction shapeFunction : shapeFunctions)
			shapeFunction.setGlobalIndex(i++);
	}
	
	private void assemblePressureFunctions(final int polynomialDegree)
	{
		
		for (final TPCell cell : getCells())
		{
			for (int i = 0; i < Math.pow(polynomialDegree + 1, getDimension()); i++)
			{
				final QkQkFunction shapeFunction = new QkQkFunction(new ContinuousTPShapeFunction(cell,
				                                                                                  polynomialDegree,
				                                                                                  i));
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (final TPCell ce : shapeFunction.getCells())
					supportOnCell.put(ce, shapeFunction);
				for (final TPFace face : shapeFunction.getFaces())
					supportOnFace.put(face, shapeFunction);
			}
		}
	}
	
	private void assembleVelocityFunctions(final int polynomialDegree)
	{
		
		for (final TPCell cell : getCells())
		{
			for (int i = 0; i < Math.pow(polynomialDegree + 2, getDimension()) * getDimension(); i++)
			{
				final QkQkFunction shapeFunction = new QkQkFunction(new ContinuousTPVectorFunction(cell,
				                                                                                   polynomialDegree +
					                                                                                   1,
				                                                                                   i));
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (final TPCell ce : shapeFunction.getCells())
					supportOnCell.put(ce, shapeFunction);
				for (final TPFace face : shapeFunction.getFaces())
					supportOnFace.put(face, shapeFunction);
			}
		}
	}
	
	public void setVelocityBoundaryValues(final VectorFunction boundaryValues)
	{
		final MixedFunction boundaryMixed = new MixedFunction(boundaryValues);
		setBoundaryValues(boundaryMixed);
	}
	
	public void setPressureBoundaryValues(final ScalarFunction boundaryValues)
	{
		final MixedFunction boundaryMixed = new MixedFunction(boundaryValues);
		setBoundaryValues(boundaryMixed);
	}
	
	public static void setVelocityBoundaryValues(final ScalarFunction indicatorFunction, final SparseMatrix mad)
	{
		throw new UnsupportedOperationException("slkdjfhl");
	}
	
	/*public void setVelocityBoundaryValues(SparseMatrix s)
	{
		setVelocityBoundaryValues(ScalarFunction.constantFunction(1), s);
	}
	
	public void setVelocityBoundaryValues(VectorFunction boundaryValues, DenseVector d)
	{
		setVelocityBoundaryValues(boundaryValues, ScalarFunction.constantFunction(1), d);
	}
	
	public void setVelocityBoundaryValues(ScalarFunction indicatorFunction, SparseMatrix s)
	{
		forEachBoundaryFace(F ->
		{
			if (TPFaceIntegral.integrateNonTensorProduct(indicatorFunction::value, F,
				QuadratureRule1D.Gauss5) > 0)
			{
				for (QkQkFunction shapeFunction : getShapeFunctionsWithSupportOnFace(F))
				{
					if (shapeFunction.hasVelocityFunction())
					{
						if (F.isOnFace(shapeFunction.getVelocityShapeFunction().getNodeFunctionalPoint()))
						{
							int shapeFunctionIndex = shapeFunction.getGlobalIndex();
							s.deleteColumn(shapeFunctionIndex);
							s.deleteRow(shapeFunctionIndex);
							s.set(1, shapeFunctionIndex, shapeFunctionIndex);
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
		forEachBoundaryFace(F ->
		{
			if (TPFaceIntegral.integrateNonTensorProduct(indicatorFunction::value,
				F, QuadratureRule1D.Gauss5) > 0)
			{
				for (QkQkFunction shapeFunction : getShapeFunctionsWithSupportOnFace(F))
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
			
		});
	}
	
	public void setPressureBoundaryValues(ScalarFunction boundaryValues)
	{
		MixedFunction boundaryMixed = new MixedFunction(boundaryValues);
		forEachBoundaryFace(face ->
		{
			for (QkQkFunction shapeFunction : getShapeFunctionsWithSupportOnFace(face))
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
			
			
		});
	}
	*/
}
