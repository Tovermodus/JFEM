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

public class TaylorHoodSpace
	extends CartesianGridSpace<QkQkFunction, MixedValue, MixedGradient, MixedHessian>
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
	
	public int getVelocitySize()
	{
		return (int) getShapeFunctions()
			.stream()
			.filter(ComposeMixedShapeFunction::hasVelocityFunction)
			.count();
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
					addFunctionToCell(shapeFunction, ce);
				for (final TPFace face : shapeFunction.getFaces())
					getFaceSupportMapping().put(face, shapeFunction);
			}
		}
	}
	
	private void assembleVelocityFunctions(final int polynomialDegree)
	{
		
		for (final TPCell cell : getCells())
		{
			for (int i = 0; i < Math.pow(polynomialDegree + 2, getDimension()) * getDimension(); i++)
			{
				final QkQkFunction shapeFunction
					= new QkQkFunction(new ContinuousTPVectorFunction(cell,
					                                                  polynomialDegree + 1,
					                                                  i));
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (final TPCell ce : shapeFunction.getCells())
					addFunctionToCell(shapeFunction, ce);
				for (final TPFace face : shapeFunction.getFaces())
					getFaceSupportMapping().put(face, shapeFunction);
			}
		}
	}
	
	public void setVelocityBoundaryValues(final VectorFunction boundaryValues)
	{
		final MixedFunction boundaryMixed = new ComposedMixedFunction(boundaryValues);
		setBoundaryValues(boundaryMixed);
	}
	
	public void setPressureBoundaryValues(final ScalarFunction boundaryValues)
	{
		final MixedFunction boundaryMixed = new ComposedMixedFunction(boundaryValues);
		setBoundaryValues(boundaryMixed);
	}
	
	public static void setVelocityBoundaryValues(final ScalarFunction indicatorFunction, final SparseMatrix mad)
	{
		throw new UnsupportedOperationException("slkdjfhl");
	}
}
