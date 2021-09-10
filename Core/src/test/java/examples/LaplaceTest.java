package examples;

import basic.*;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import org.junit.Test;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class LaplaceTest
{
	@Test(timeout = 20000)
	public void testDGConvergence()
	{
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int mult = 4;
		final List<Integer> cells = List.of(mult * 2, mult * 3);
		
		final double penalty = 1000;
		
		final TPCellIntegral<TPShapeFunction> gg =
			new TPCellIntegralViaReferenceCell<>(TPCellIntegral.GRAD_GRAD);
		final TPFaceIntegral<TPShapeFunction> jj =
			new TPFaceIntegralViaReferenceFace<>(penalty, TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
		final TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
			                              TPRightHandSideIntegral.VALUE);
		final TPBoundaryFaceIntegral<TPShapeFunction> boundaryValues =
			new TPBoundaryFaceIntegral<>(LaplaceReferenceSolution.scalarBoundaryValues(penalty),
			                             TPBoundaryFaceIntegral.VALUE);
		
		final ArrayList<CellIntegral<TPCell, TPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		final ArrayList<FaceIntegral<TPFace, TPShapeFunction>> faceIntegrals =
			new ArrayList<>();
		final ArrayList<RightHandSideIntegral<TPCell, TPShapeFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		final ArrayList<BoundaryRightHandSideIntegral<TPFace, TPShapeFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		boundaryFaceIntegrals.add(boundaryValues);
		rightHandSideIntegrals.add(rightHandSideIntegral);
		cellIntegrals.add(gg);
		faceIntegrals.add(jj);
		
		final int polynomialDegree = 2;
		final TPFESpace grid = new TPFESpace(start, end, cells);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = false;
		final Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		final ScalarFESpaceFunction<TPShapeFunction> solut =
			new ScalarFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		assertTrue(ConvergenceOrderEstimator.normL2Difference(solut,
		                                                      LaplaceReferenceSolution.scalarReferenceSolution(),
		                                                      grid.generatePlotPoints(20)) < 1e-2);
	}
	
	@Test(timeout = 20000)
	public void testContinuousConvergence()
	{
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int mult = 4;
		final List<Integer> cells = List.of(mult * 2, mult * 3);
		
		final TPCellIntegral<ContinuousTPShapeFunction> gg =
			new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		final TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
			                              TPRightHandSideIntegral.VALUE);
		
		final ArrayList<CellIntegral<TPCell, ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		final ArrayList<FaceIntegral<TPFace, ContinuousTPShapeFunction>> faceIntegrals =
			new ArrayList<>();
		final ArrayList<RightHandSideIntegral<TPCell, ContinuousTPShapeFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		final ArrayList<BoundaryRightHandSideIntegral<TPFace, ContinuousTPShapeFunction>> boundaryFaceIntegrals
			=
			new ArrayList<>();
		
		rightHandSideIntegrals.add(rightHandSideIntegral);
		cellIntegrals.add(gg);
		
		final int polynomialDegree = 2;
		final ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end, cells);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setBoundaryValues(LaplaceReferenceSolution.scalarBoundaryValues());
		
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = false;
		final Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> solut =
			new ScalarFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		assertTrue(ConvergenceOrderEstimator.normL2Difference(solut,
		                                                      LaplaceReferenceSolution.scalarReferenceSolution(),
		                                                      grid.generatePlotPoints(20)) < 1e-2);
	}
}
