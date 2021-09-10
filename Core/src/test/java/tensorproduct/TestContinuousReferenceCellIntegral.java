package tensorproduct;

import basic.CellIntegral;
import basic.PerformanceArguments;
import basic.RightHandSideIntegral;
import basic.ScalarFunction;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.Matrix;
import org.junit.Test;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.Cell1D;
import tensorproduct.geometry.TPCell;

import java.util.ArrayList;

import static org.junit.Assert.assertTrue;

public class TestContinuousReferenceCellIntegral
{
	@Test
	public void test2DGradGradMatrix()
	{
		final Matrix referenceMatrix;
		final CoordinateVector start = CoordinateVector.fromValues(-1, -10);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 3;
		
		ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end,
		                                                   Ints.asList(3, 5));
		TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
		                                                                    TPCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell, ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1),
			                              TPRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell, ContinuousTPShapeFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new ContinuousTPFESpace(start, end,
		                               Ints.asList(3, 5));
		gg = new TPCellIntegralViaReferenceCell<>(1,
		                                          TPCellIntegral.GRAD_GRAD);
		cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1),
			                              TPRightHandSideIntegral.VALUE);
		rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix()));
	}
	
	@Test
	public void test3DGradGradMatrix()
	{
		final Matrix referenceMatrix;
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(7, 4, 119);
		final int polynomialDegree = 2;
		
		ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end,
		                                                   Ints.asList(4, 7, 2));
		TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
		                                                                    TPCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell, ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1),
			                              TPRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell, ContinuousTPShapeFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new ContinuousTPFESpace(start, end,
		                               Ints.asList(4, 7, 2));
		gg = new TPCellIntegralViaReferenceCell<>(1,
		                                          TPCellIntegral.GRAD_GRAD);
		cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1),
			                              TPRightHandSideIntegral.VALUE);
		rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix()));
	}
	
	@Test
	public void test2DGradGrad()
	{
		final CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(3, 0),
		                                          CoordinateVector.fromValues(5, 4), new IntCoordinates(1, 1));
		final TPCell c = g.cells.get(0);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					final ContinuousTPShapeFunction f1 = new ContinuousTPShapeFunction(c, k, i);
					final ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c, k, j);
					final TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(
						ScalarFunction.constantFunction(1),
						TPCellIntegral.GRAD_GRAD);
					final TPCellIntegral<ContinuousTPShapeFunction> gg2
						= new TPCellIntegralViaReferenceCell<>(1,
						                                       TPCellIntegral.GRAD_GRAD);
					assertTrue(Math.abs(
						gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						                                                              f1,
						                                                              f2)) < 1e-12);
				}
	}
	
	@Test
	public void test3DGradGrad()
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		final ArrayList<Cell1D> cell1DList = new ArrayList<>();
		final CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0, 0, 0),
		                                          CoordinateVector.fromValues(2, 1, 1), new IntCoordinates(1, 1, 1));
		final TPCell c = g.cells.get(0);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					final ContinuousTPShapeFunction f1 = new ContinuousTPShapeFunction(c, k, i);
					final ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c, k, j);
					final TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(
						ScalarFunction.constantFunction(1),
						TPCellIntegral.GRAD_GRAD);
					final TPCellIntegral<ContinuousTPShapeFunction> gg2
						= new TPCellIntegralViaReferenceCell<>(1,
						                                       TPCellIntegral.GRAD_GRAD);
					assertTrue(Math.abs(
						gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						                                                              f1,
						                                                              f2)) < 1e-12);
				}
	}
	
	@Test
	public void test3DValueValue()
	{
		final CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0, 0, 0),
		                                          CoordinateVector.fromValues(2, 1, 1), new IntCoordinates(1, 1, 1));
		final TPCell c = g.cells.get(0);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					final ContinuousTPShapeFunction f1 = new ContinuousTPShapeFunction(c, k, i);
					final ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c, k, j);
					final TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(
						ScalarFunction.constantFunction(1),
						TPCellIntegral.VALUE_VALUE);
					final TPCellIntegral<ContinuousTPShapeFunction> gg2
						= new TPCellIntegralViaReferenceCell<>(1,
						                                       TPCellIntegral.VALUE_VALUE);
					assertTrue(Math.abs(
						gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						                                                              f1,
						                                                              f2)) < 1e-12);
				}
	}
}
