package distorted;

import basic.DoubleCompare;
import distorted.geometry.DistortedCell;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.junit.jupiter.api.Test;
import tensorproduct.TPCellIntegral;
import tensorproduct.TPShapeFunction;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.TPCell;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class DistortedCellIntegralTest
{
	private static DistortedCell getNonDistortedCell2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1);
		vertices[3] = CoordinateVector.fromValues(0, 1);
		return new DistortedCell(vertices);
	}
	
	private static DistortedCell getScaled2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(4, 0);
		vertices[2] = CoordinateVector.fromValues(4, 2);
		vertices[3] = CoordinateVector.fromValues(0, 2);
		return new DistortedCell(vertices);
	}
	
	private static DistortedCell getRotated2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(1, -1);
		vertices[2] = CoordinateVector.fromValues(2, 0);
		vertices[3] = CoordinateVector.fromValues(1, 1);
		return new DistortedCell(vertices);
	}
	
	@Test
	public void testNonDistortedValue()
	{
		final DistortedCell cell = getNonDistortedCell2D();
		final CartesianGrid grid = new CartesianGrid(CoordinateVector.fromValues(0, 0),
		                                             CoordinateVector.fromValues(1, 1),
		                                             new IntCoordinates(1, 1));
		final TPCell cellRef = grid.cells.get(0);
		final DistortedCellIntegral gradGrad = new DistortedCellIntegral(DistortedCellIntegral.VALUE_VALUE);
		final TPCellIntegral<TPShapeFunction> gradGradRef = new TPCellIntegral<>(TPCellIntegral.VALUE_VALUE);
		final List<TPShapeFunction> shapeFunctionsRef = new ArrayList<>();
		final List<DistortedShapeFunction> shapeFunctions = new ArrayList<>();
		final int polynomialDegree = 1;
		for (int i = 0; i < Math.pow(polynomialDegree + 1, 2); i++)
		{
			shapeFunctionsRef.add(new TPShapeFunction(cellRef, polynomialDegree, i));
			shapeFunctions.add(new DistortedShapeFunction(cell, polynomialDegree, i));
		}
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				assertTrue(DoubleCompare.almostEqual(
					gradGradRef.evaluateCellIntegral(cellRef, shapeFunctionsRef.get(i),
					                                 shapeFunctionsRef.get(j))
					, gradGrad.evaluateCellIntegral(cell, shapeFunctions.get(i),
					                                shapeFunctions.get(j))));
			}
		}
	}
	
	@Test
	public void testScaledValue()
	{
		final DistortedCell cell = getScaled2D();
		final CartesianGrid grid = new CartesianGrid(CoordinateVector.fromValues(0, 0),
		                                             CoordinateVector.fromValues(4, 2),
		                                             new IntCoordinates(1, 1));
		final TPCell cellRef = grid.cells.get(0);
		final DistortedCellIntegral gradGrad = new DistortedCellIntegral(DistortedCellIntegral.VALUE_VALUE);
		final TPCellIntegral<TPShapeFunction> gradGradRef = new TPCellIntegral<>(TPCellIntegral.VALUE_VALUE);
		final List<TPShapeFunction> shapeFunctionsRef = new ArrayList<>();
		final List<DistortedShapeFunction> shapeFunctions = new ArrayList<>();
		final int polynomialDegree = 1;
		for (int i = 0; i < Math.pow(polynomialDegree + 1, 2); i++)
		{
			shapeFunctionsRef.add(new TPShapeFunction(cellRef, polynomialDegree, i));
			shapeFunctions.add(new DistortedShapeFunction(cell, polynomialDegree, i));
		}
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				assertTrue(DoubleCompare.almostEqual(
					gradGradRef.evaluateCellIntegral(cellRef, shapeFunctionsRef.get(i),
					                                 shapeFunctionsRef.get(j))
					, gradGrad.evaluateCellIntegral(cell, shapeFunctions.get(i),
					                                shapeFunctions.get(j))));
			}
		}
	}
	
	@Test
	public void testRotatedValue()
	{
		final DistortedCell cell = getRotated2D();
		final CartesianGrid grid = new CartesianGrid(CoordinateVector.fromValues(0, 0),
		                                             CoordinateVector.fromValues(Math.sqrt(2), Math.sqrt(2)),
		                                             new IntCoordinates(1, 1));
		final TPCell cellRef = grid.cells.get(0);
		final DistortedCellIntegral gradGrad = new DistortedCellIntegral(DistortedCellIntegral.VALUE_VALUE);
		final TPCellIntegral<TPShapeFunction> gradGradRef = new TPCellIntegral<>(TPCellIntegral.VALUE_VALUE);
		final List<TPShapeFunction> shapeFunctionsRef = new ArrayList<>();
		final List<DistortedShapeFunction> shapeFunctions = new ArrayList<>();
		final int polynomialDegree = 1;
		for (int i = 0; i < Math.pow(polynomialDegree + 1, 2); i++)
		{
			shapeFunctionsRef.add(new TPShapeFunction(cellRef, polynomialDegree, i));
			shapeFunctions.add(new DistortedShapeFunction(cell, polynomialDegree, i));
		}
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				assertTrue(DoubleCompare.almostEqual(
					gradGradRef.evaluateCellIntegral(cellRef, shapeFunctionsRef.get(i),
					                                 shapeFunctionsRef.get(j))
					, gradGrad.evaluateCellIntegral(cell, shapeFunctions.get(i),
					                                shapeFunctions.get(j))));
			}
		}
	}
	
	@Test
	public void testNonDistortedGradient()
	{
		final DistortedCell cell = getNonDistortedCell2D();
		final CartesianGrid grid = new CartesianGrid(CoordinateVector.fromValues(0, 0),
		                                             CoordinateVector.fromValues(1, 1),
		                                             new IntCoordinates(1, 1));
		final TPCell cellRef = grid.cells.get(0);
		final DistortedCellIntegral gradGrad = new DistortedCellIntegral(DistortedCellIntegral.GRAD_GRAD);
		final TPCellIntegral<TPShapeFunction> gradGradRef = new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		final List<TPShapeFunction> shapeFunctionsRef = new ArrayList<>();
		final List<DistortedShapeFunction> shapeFunctions = new ArrayList<>();
		final int polynomialDegree = 1;
		for (int i = 0; i < Math.pow(polynomialDegree + 1, 2); i++)
		{
			shapeFunctionsRef.add(new TPShapeFunction(cellRef, polynomialDegree, i));
			shapeFunctions.add(new DistortedShapeFunction(cell, polynomialDegree, i));
		}
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				assertTrue(DoubleCompare.almostEqual(
					gradGradRef.evaluateCellIntegral(cellRef, shapeFunctionsRef.get(i),
					                                 shapeFunctionsRef.get(j))
					, gradGrad.evaluateCellIntegral(cell, shapeFunctions.get(i),
					                                shapeFunctions.get(j))));
			}
		}
	}
	
	@Test
	public void testScaledGradient()
	{
		final DistortedCell cell = getScaled2D();
		final CartesianGrid grid = new CartesianGrid(CoordinateVector.fromValues(0, 0),
		                                             CoordinateVector.fromValues(4, 2),
		                                             new IntCoordinates(1, 1));
		final TPCell cellRef = grid.cells.get(0);
		final DistortedCellIntegral gradGrad = new DistortedCellIntegral(DistortedCellIntegral.GRAD_GRAD);
		final TPCellIntegral<TPShapeFunction> gradGradRef = new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		final List<TPShapeFunction> shapeFunctionsRef = new ArrayList<>();
		final List<DistortedShapeFunction> shapeFunctions = new ArrayList<>();
		final int polynomialDegree = 1;
		for (int i = 0; i < Math.pow(polynomialDegree + 1, 2); i++)
		{
			shapeFunctionsRef.add(new TPShapeFunction(cellRef, polynomialDegree, i));
			shapeFunctions.add(new DistortedShapeFunction(cell, polynomialDegree, i));
		}
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				assertTrue(DoubleCompare.almostEqual(
					gradGradRef.evaluateCellIntegral(cellRef, shapeFunctionsRef.get(i),
					                                 shapeFunctionsRef.get(j))
					, gradGrad.evaluateCellIntegral(cell, shapeFunctions.get(i),
					                                shapeFunctions.get(j))));
			}
		}
	}
	
	@Test
	public void testRotatedGradient()
	{
		final DistortedCell cell = getRotated2D();
		final CartesianGrid grid = new CartesianGrid(CoordinateVector.fromValues(0, 0),
		                                             CoordinateVector.fromValues(Math.sqrt(2), Math.sqrt(2)),
		                                             new IntCoordinates(1, 1));
		final TPCell cellRef = grid.cells.get(0);
		final DistortedCellIntegral gradGrad = new DistortedCellIntegral(DistortedCellIntegral.GRAD_GRAD);
		final TPCellIntegral<TPShapeFunction> gradGradRef = new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		final List<TPShapeFunction> shapeFunctionsRef = new ArrayList<>();
		final List<DistortedShapeFunction> shapeFunctions = new ArrayList<>();
		final int polynomialDegree = 1;
		for (int i = 0; i < Math.pow(polynomialDegree + 1, 2); i++)
		{
			shapeFunctionsRef.add(new TPShapeFunction(cellRef, polynomialDegree, i));
			shapeFunctions.add(new DistortedShapeFunction(cell, polynomialDegree, i));
		}
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				assertTrue(DoubleCompare.almostEqual(
					gradGradRef.evaluateCellIntegral(cellRef, shapeFunctionsRef.get(i),
					                                 shapeFunctionsRef.get(j))
					, gradGrad.evaluateCellIntegral(cell, shapeFunctions.get(i),
					                                shapeFunctions.get(j))));
			}
		}
	}
}
