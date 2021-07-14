package tensorproduct;

import basic.CellIntegral;
import basic.PerformanceArguments;
import basic.RightHandSideIntegral;
import basic.ScalarFunction;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.Matrix;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class TestContinuousReferenceCellIntegral
{
	@Test
	public void test2DGradGradMatrix()
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		long startTime = System.nanoTime();
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -10);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 3;
		
		ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end,
			Ints.asList(3, 5), polynomialDegree);
		TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD,
			false);
		ArrayList<CellIntegral<TPCell, ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1), TPRightHandSideIntegral.VALUE,
				true);
		ArrayList<RightHandSideIntegral<TPCell, ContinuousTPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		System.out.println(grid.getCells());
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new ContinuousTPFESpace(start, end,
			Ints.asList(3, 5), polynomialDegree);
		gg = new TPCellIntegralViaReferenceCell<>(1,
			TPCellIntegral.GRAD_GRAD,
			false);
		cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1), TPRightHandSideIntegral.VALUE,
				true);
		rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix(), 1e-12));
		
	}
	
	@Test
	public void test3DGradGradMatrix()
	{
		long startTime = System.nanoTime();
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -1, -1);
		CoordinateVector end = CoordinateVector.fromValues(7, 4, 119);
		int polynomialDegree = 2;
		
		ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end,
			Ints.asList(4, 7, 2), polynomialDegree);
		TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD,
			false);
		ArrayList<CellIntegral<TPCell, ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1), TPRightHandSideIntegral.VALUE,
				true);
		ArrayList<RightHandSideIntegral<TPCell, ContinuousTPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new ContinuousTPFESpace(start, end,
			Ints.asList(4, 7, 2), polynomialDegree);
		gg = new TPCellIntegralViaReferenceCell<>(1,
			TPCellIntegral.GRAD_GRAD,
			false);
		cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1), TPRightHandSideIntegral.VALUE,
				true);
		rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix(), 1e-12));
		
	}
	
	@Test
	public void test2DGradGrad()
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		ArrayList<Cell1D> cell1DList = new ArrayList<>();
		cell1DList.add(new Cell1D(3, 5)); //2/3, 5/6, 10/9, 17/12, 26/15, 37/18
		cell1DList.add(new Cell1D(0, 4));// 2/3, 13/18 5/6 29/30 10/9 53/42
		TPCell c = new TPCell(cell1DList);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					ContinuousTPShapeFunction f1 = new ContinuousTPShapeFunction(c, k, i);
					ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c, k, j);
					TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
						TPCellIntegral.GRAD_GRAD,
						false);
					TPCellIntegral<ContinuousTPShapeFunction> gg2 = new TPCellIntegralViaReferenceCell<>(1,
						TPCellIntegral.GRAD_GRAD, false);
					assertTrue(Math.abs(gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						f1, f2)) < 1e-12, "Difference of " + Math.abs(gg.evaluateCellIntegral(c, f1,
						f2) - gg2.evaluateCellIntegral(c,
						f1, f2)) + " for polynomial degree " + k + " and functions " + i + " and " + j);
				}
	}
	
	@Test
	public void test3DGradGrad()
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		ArrayList<Cell1D> cell1DList = new ArrayList<>();
		cell1DList.add(new Cell1D(0, 2));
		cell1DList.add(new Cell1D(0, 1));
		cell1DList.add(new Cell1D(0, 1));
		TPCell c = new TPCell(cell1DList);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					System.out.println("--------------" + i + " " + j + " k " + k);
					ContinuousTPShapeFunction f1 = new ContinuousTPShapeFunction(c, k, i);
					ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c, k, j);
					System.out.println(f1);
					System.out.println(f2);
					TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
						TPCellIntegral.GRAD_GRAD,
						false);
					TPCellIntegral<ContinuousTPShapeFunction> gg2 = new TPCellIntegralViaReferenceCell<>(1,
						TPCellIntegral.GRAD_GRAD, false);
					System.out.println("normal int " + gg.evaluateCellIntegral(c, f1, f2));
					System.out.println("reference int " + gg2.evaluateCellIntegral(c, f1, f2));
					assertTrue(Math.abs(gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						f1, f2)) < 1e-12, "Difference of " + Math.abs(gg.evaluateCellIntegral(c, f1,
						f2) - gg2.evaluateCellIntegral(c,
						f1, f2)) + " for polynomial degree " + k + " and functions " + i + " and " + j);
				}
	}
	
	@Test
	public void test3DValueValue()
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		ArrayList<Cell1D> cell1DList = new ArrayList<>();
		cell1DList.add(new Cell1D(0, 2));
		cell1DList.add(new Cell1D(0, 1));
		cell1DList.add(new Cell1D(0, 1));
		TPCell c = new TPCell(cell1DList);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					System.out.println("--------------" + i + " " + j + " k " + k);
					ContinuousTPShapeFunction f1 = new ContinuousTPShapeFunction(c, k, i);
					ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c, k, j);
					System.out.println(f1);
					System.out.println(f2);
					TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
						TPCellIntegral.VALUE_VALUE,
						false);
					TPCellIntegral<ContinuousTPShapeFunction> gg2 = new TPCellIntegralViaReferenceCell<>(1,
						TPCellIntegral.VALUE_VALUE, false);
					assertTrue(Math.abs(gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						f1, f2)) < 1e-12, "Difference of " + Math.abs(gg.evaluateCellIntegral(c, f1,
						f2) - gg2.evaluateCellIntegral(c,
						f1, f2)) + " for polynomial degree " + k + " and functions " + i + " and " + j);
				}
	}
}
