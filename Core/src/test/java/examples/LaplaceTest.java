package examples;

import basic.*;
import linalg.*;
import org.junit.jupiter.api.Test;
import systems.SystemMixedCellIntegral;
import systems.SystemParameters;
import systems.SystemShapeFunction;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class LaplaceTest
{
	@Test
	public void testDGConvergence()
	{
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int mult = 4;
		List<Integer> cells = List.of(mult*2,mult*3);
		
		double penalty = 1000;
		
		TPCellIntegral<TPShapeFunction> gg =
			new TPCellIntegralViaReferenceCell<>(TPCellIntegral.GRAD_GRAD);
		TPFaceIntegral<TPShapeFunction> jj =
			new TPFaceIntegralViaReferenceFace<>(penalty,TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
		TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
				TPRightHandSideIntegral.VALUE);
		TPBoundaryFaceIntegral<TPShapeFunction> boundaryValues =
			new TPBoundaryFaceIntegral<>(LaplaceReferenceSolution.scalarBoundaryValues(penalty),
				TPBoundaryFaceIntegral.VALUE);
		
		ArrayList<CellIntegral<TPCell,TPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		ArrayList<FaceIntegral< TPFace,TPShapeFunction>> faceIntegrals =
			new ArrayList<>();
		ArrayList<RightHandSideIntegral<TPCell, TPShapeFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		ArrayList<BoundaryRightHandSideIntegral<TPFace, TPShapeFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		boundaryFaceIntegrals.add(boundaryValues);
		rightHandSideIntegrals.add(rightHandSideIntegral);
		cellIntegrals.add(gg);
		faceIntegrals.add(jj);
		
		int polynomialDegree = 2;
		TPFESpace grid = new TPFESpace(start, end, cells, polynomialDegree);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		
		IterativeSolver it = new IterativeSolver();
		it.showProgress = false;
		Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		ScalarFESpaceFunction<TPShapeFunction> solut =
			new ScalarFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		assertTrue(ConvergenceOrderEstimator.normL2Difference(solut,
			LaplaceReferenceSolution.scalarReferenceSolution(), grid.generatePlotPoints(20)) < 1e-2);
	}
	
	@Test
	public void testContinuousConvergence()
	{
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int mult = 4;
		List<Integer> cells = List.of(mult*2,mult*3);
		
		
		TPCellIntegral<ContinuousTPShapeFunction> gg =
			new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
				TPRightHandSideIntegral.VALUE);
		
		ArrayList<CellIntegral<TPCell,ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		ArrayList<FaceIntegral< TPFace,ContinuousTPShapeFunction>> faceIntegrals =
			new ArrayList<>();
		ArrayList<RightHandSideIntegral<TPCell, ContinuousTPShapeFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		ArrayList<BoundaryRightHandSideIntegral<TPFace, ContinuousTPShapeFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		rightHandSideIntegrals.add(rightHandSideIntegral);
		cellIntegrals.add(gg);
		
		int polynomialDegree = 2;
		ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end, cells, polynomialDegree);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setBoundaryValues(LaplaceReferenceSolution.scalarBoundaryValues());
		
		IterativeSolver it = new IterativeSolver();
		it.showProgress = false;
		Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		ScalarFESpaceFunction<ContinuousTPShapeFunction> solut =
			new ScalarFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		assertTrue(ConvergenceOrderEstimator.normL2Difference(solut,
			LaplaceReferenceSolution.scalarReferenceSolution(), grid.generatePlotPoints(20)) < 1e-2);
	}
}
