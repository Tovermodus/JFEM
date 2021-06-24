package tensorproduct;

import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.Matrix;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class TestReferenceIntegral
{
	@Test
	public void test2DGradGradMatrix()
	{
		PerformanceArguments.createInstance(true,12,true);
		
		long startTime = System.nanoTime();
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1,-10);
		CoordinateVector end = CoordinateVector.fromValues(1,1);
		int polynomialDegree = 3;
		
		TPFESpace grid = new TPFESpace(start,end,
			Ints.asList(6,6),polynomialDegree);
		TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD,
			false);
		ArrayList<CellIntegral<TPCell,TPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1),TPRightHandSideIntegral.VALUE,
				true);
		ArrayList<RightHandSideIntegral<TPCell,TPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new TPFESpace(start,end,
			Ints.asList(6,6,6),polynomialDegree);
		gg = new TPCellIntegralViaReferenceCell<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD,
			false);
		cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1),TPRightHandSideIntegral.VALUE,
				true);
		rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix()));
		
	}
	@Test
	public void test3DGradGradMatrix()
	{
		PerformanceArguments.createInstance(true,12,true);
		
		long startTime = System.nanoTime();
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1,-1,-1);
		CoordinateVector end = CoordinateVector.fromValues(7,4,1);
		int polynomialDegree = 1;
		
		TPFESpace grid = new TPFESpace(start,end,
			Ints.asList(4,4,4),polynomialDegree);
		TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD,
			false);
		ArrayList<CellIntegral<TPCell,TPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1),TPRightHandSideIntegral.VALUE,
				true);
		ArrayList<RightHandSideIntegral<TPCell,TPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new TPFESpace(start,end,
			Ints.asList(4,4,4),polynomialDegree);
		gg = new TPCellIntegralViaReferenceCell<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD,
			false);
		cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1),TPRightHandSideIntegral.VALUE,
				true);
		rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix()));
		
	}
	public void test3DGradGrad()
	{
	
	}
}
