package tensorproduct;

import basic.CellIntegral;
import basic.PerformanceArguments;
import basic.RightHandSideIntegral;
import basic.ScalarFunction;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.Matrix;
import org.junit.jupiter.api.Test;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.Cell1D;
import tensorproduct.geometry.TPCell;

import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class TestContinuousReferenceCellIntegral
{
	@Test
	public void test2DGradGradMatrix()
	{
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -10);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 3;
		
		ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end,
			Ints.asList(3, 5));
		TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell, ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1), TPRightHandSideIntegral.VALUE);
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
			Ints.asList(3, 5));
		gg = new TPCellIntegralViaReferenceCell<>(1,
			TPCellIntegral.GRAD_GRAD);
		cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1), TPRightHandSideIntegral.VALUE);
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
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -1, -1);
		CoordinateVector end = CoordinateVector.fromValues(7, 4, 119);
		int polynomialDegree = 2;
		
		ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end,
			Ints.asList(4, 7, 2));
		TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell, ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1), TPRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell, ContinuousTPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
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
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1), TPRightHandSideIntegral.VALUE);
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
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(3,0),
			CoordinateVector.fromValues(5,4),new IntCoordinates(1,1));
		TPCell c = g.cells.get(0);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					ContinuousTPShapeFunction f1 = new ContinuousTPShapeFunction(c, k, i);
					ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c, k, j);
					TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
						TPCellIntegral.GRAD_GRAD);
					TPCellIntegral<ContinuousTPShapeFunction> gg2 = new TPCellIntegralViaReferenceCell<>(1,
						TPCellIntegral.GRAD_GRAD);
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
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0,0,0),
			CoordinateVector.fromValues(2,1,1),new IntCoordinates(1,1,1));
		TPCell c = g.cells.get(0);
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
						TPCellIntegral.GRAD_GRAD);
					TPCellIntegral<ContinuousTPShapeFunction> gg2 = new TPCellIntegralViaReferenceCell<>(1,
						TPCellIntegral.GRAD_GRAD);
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
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0,0,0),
			CoordinateVector.fromValues(2,1,1),new IntCoordinates(1,1,1));
		TPCell c = g.cells.get(0);
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
						TPCellIntegral.VALUE_VALUE);
					TPCellIntegral<ContinuousTPShapeFunction> gg2 = new TPCellIntegralViaReferenceCell<>(1,
						TPCellIntegral.VALUE_VALUE);
					assertTrue(Math.abs(gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						f1, f2)) < 1e-12, "Difference of " + Math.abs(gg.evaluateCellIntegral(c, f1,
						f2) - gg2.evaluateCellIntegral(c,
						f1, f2)) + " for polynomial degree " + k + " and functions " + i + " and " + j);
				}
	}
}
