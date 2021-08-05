package tensorproduct;

import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.Matrix;
import org.junit.jupiter.api.Test;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.TPCell;

import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class TestReferenceCellIntegral
{
	@Test
	public void test2DGradGradMatrix()
	{
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -10);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 3;
		
		TPFESpace grid = new TPFESpace(start, end,
			Ints.asList(3, 5));
		TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell, TPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1), TPRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell, TPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new TPFESpace(start, end,
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
		
		TPFESpace grid = new TPFESpace(start, end,
			Ints.asList(4, 7, 2));
		TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell, TPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1), TPRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell, TPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new TPFESpace(start, end,
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
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(3,0), CoordinateVector.fromValues(5,
			4),new IntCoordinates(1,1));
		TPCell c = g.cells.get(0);
		for(int k = 0; k < 5; k++)
		for (int i = 0; i < (k + 1) * (k + 1); i++)
			for (int j = 0; j < (k + 1) * (k + 1); j++)
			{
				TPShapeFunction f1 = new TPShapeFunction(c, k, i);
				TPShapeFunction f2 = new TPShapeFunction(c, k, j);
				TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
					TPCellIntegral.GRAD_GRAD);
				TPCellIntegral<TPShapeFunction> gg2 = new TPCellIntegralViaReferenceCell<>(1,
					TPCellIntegral.GRAD_GRAD);
				assertTrue(Math.abs(gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
					f1, f2)) < 1e-12, "Difference of "+Math.abs(gg.evaluateCellIntegral(c, f1,
					f2) - gg2.evaluateCellIntegral(c,
					f1, f2))+ " for polynomial degree " + k + " and functions " + i+" and " + j);
			}
	}
	
	@Test
	public void test3DGradGrad()
	{
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0,0,0),
			CoordinateVector.fromValues(2,1,1),new IntCoordinates(1,1,1));
		TPCell c = g.cells.get(0);
		for(int k = 0; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					TPShapeFunction f1 = new TPShapeFunction(c, k, i);
					TPShapeFunction f2 = new TPShapeFunction(c, k, j);
					TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
						TPCellIntegral.GRAD_GRAD);
					TPCellIntegral<TPShapeFunction> gg2 = new TPCellIntegralViaReferenceCell<>(1,
						TPCellIntegral.GRAD_GRAD);
					assertTrue(Math.abs(gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						f1, f2)) < 1e-12, "Difference of "+Math.abs(gg.evaluateCellIntegral(c, f1,
						f2) - gg2.evaluateCellIntegral(c,
						f1, f2))+ " for polynomial degree " + k + " and functions " + i+" and " + j);
				}
	}
	@Test
	public void test3DValueValue()
	{
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0,0,0),
			CoordinateVector.fromValues(2,1,1),new IntCoordinates(1,1,1));
		TPCell c = g.cells.get(0);
		for(int k = 0; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					TPShapeFunction f1 = new TPShapeFunction(c, k, i);
					TPShapeFunction f2 = new TPShapeFunction(c, k, j);
					TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
						TPCellIntegral.VALUE_VALUE);
					TPCellIntegral<TPShapeFunction> gg2 = new TPCellIntegralViaReferenceCell<>(1,
						TPCellIntegral.VALUE_VALUE);
					assertTrue(Math.abs(gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						f1, f2)) < 1e-12, "Difference of "+Math.abs(gg.evaluateCellIntegral(c, f1,
						f2) - gg2.evaluateCellIntegral(c,
						f1, f2))+ " for polynomial degree " + k + " and functions " + i+" and " + j);
				}
	}
}
