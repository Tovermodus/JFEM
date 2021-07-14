package examples;

import basic.*;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.SparseMatrix;
import linalg.Vector;
import org.junit.jupiter.api.Test;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class VectorLaplaceTest
{
	@Test
	public void testDGConvergence()
	{
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int mult = 4;
		List<Integer> cells = List.of(mult*2,mult*3);
		
		double penalty = 4000;
		
		TPVectorCellIntegral<TPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		TPVectorFaceIntegral<TPVectorFunction> jj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(penalty),
				TPVectorFaceIntegral.VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE);
		
		TPVectorRightHandSideIntegral<TPVectorFunction> rightHandSideIntegral =
			new TPVectorRightHandSideIntegral<>(LaplaceReferenceSolution.vectorRightHandSide(),
				TPVectorRightHandSideIntegral.VALUE);
		TPVectorBoundaryFaceIntegral<TPVectorFunction> boundaryValues =
			new TPVectorBoundaryFaceIntegral<>(LaplaceReferenceSolution.vectorBoundaryValues(penalty/4),
				TPVectorBoundaryFaceIntegral.VALUE);
		
		ArrayList<CellIntegral<TPCell,TPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		ArrayList<FaceIntegral< TPFace,TPVectorFunction>> faceIntegrals =
			new ArrayList<>();
		ArrayList<RightHandSideIntegral<TPCell, TPVectorFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		ArrayList<BoundaryRightHandSideIntegral<TPFace, TPVectorFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		boundaryFaceIntegrals.add(boundaryValues);
		rightHandSideIntegrals.add(rightHandSideIntegral);
		cellIntegrals.add(gg);
		faceIntegrals.add(jj);
		
		int polynomialDegree = 2;
		TPVectorFESpace grid = new TPVectorFESpace(start, end, cells, polynomialDegree);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		
		IterativeSolver<SparseMatrix> it = new IterativeSolver<>();
		it.showProgress = false;
		Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		VectorFESpaceFunction<TPVectorFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(solut,
			LaplaceReferenceSolution.vectorReferenceSolution(), grid.generatePlotPoints(20)) < 1e-2);
	}
	
	@Test
	public void testContinuousConvergence()
	{
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int mult = 4;
		List<Integer> cells = List.of(mult*2,mult*3);
		
		double penalty = 100;
		
		TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		
		TPVectorRightHandSideIntegral<ContinuousTPVectorFunction> rightHandSideIntegral =
			new TPVectorRightHandSideIntegral<>(LaplaceReferenceSolution.vectorRightHandSide(),
				TPVectorRightHandSideIntegral.VALUE);
		
		ArrayList<CellIntegral<TPCell,ContinuousTPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		ArrayList<FaceIntegral< TPFace,ContinuousTPVectorFunction>> faceIntegrals =
			new ArrayList<>();
		ArrayList<RightHandSideIntegral<TPCell, ContinuousTPVectorFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		ArrayList<BoundaryRightHandSideIntegral<TPFace, ContinuousTPVectorFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		rightHandSideIntegrals.add(rightHandSideIntegral);
		cellIntegrals.add(gg);
		
		int polynomialDegree = 2;
		ContinuousTPFEVectorSpace grid = new ContinuousTPFEVectorSpace(start, end, cells, polynomialDegree);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setBoundaryValues(LaplaceReferenceSolution.vectorBoundaryValues());
		IterativeSolver<SparseMatrix> it = new IterativeSolver<>();
		it.showProgress = false;
		Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		VectorFESpaceFunction<ContinuousTPVectorFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(solut,
			LaplaceReferenceSolution.vectorReferenceSolution(), grid.generatePlotPoints(20)) < 1e-2);
	}
}
