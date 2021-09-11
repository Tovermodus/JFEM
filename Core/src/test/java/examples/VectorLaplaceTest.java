package examples;

import basic.*;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.Timeout;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class VectorLaplaceTest
{
	@Rule
	public Timeout testTimeout = new Timeout(20000);
	
	@Test
	public void testDGConvergence()
	{
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int mult = 4;
		final List<Integer> cells = List.of(mult * 2, mult * 3);
		
		final double penalty = 4000;
		
		final TPVectorCellIntegral<TPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		final TPVectorFaceIntegral<TPVectorFunction> jj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(penalty),
			                           TPVectorFaceIntegral.VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE);
		
		final TPVectorRightHandSideIntegral<TPVectorFunction> rightHandSideIntegral =
			new TPVectorRightHandSideIntegral<>(LaplaceReferenceSolution.vectorRightHandSide(),
			                                    TPVectorRightHandSideIntegral.VALUE);
		final TPVectorBoundaryFaceIntegral<TPVectorFunction> boundaryValues =
			new TPVectorBoundaryFaceIntegral<>(LaplaceReferenceSolution.vectorBoundaryValues(penalty / 4),
			                                   TPVectorBoundaryFaceIntegral.VALUE);
		
		final ArrayList<CellIntegral<TPCell, TPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		final ArrayList<FaceIntegral<TPFace, TPVectorFunction>> faceIntegrals =
			new ArrayList<>();
		final ArrayList<RightHandSideIntegral<TPCell, TPVectorFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		final ArrayList<BoundaryRightHandSideIntegral<TPFace, TPVectorFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		boundaryFaceIntegrals.add(boundaryValues);
		rightHandSideIntegrals.add(rightHandSideIntegral);
		cellIntegrals.add(gg);
		faceIntegrals.add(jj);
		
		final int polynomialDegree = 2;
		final TPVectorFESpace grid = new TPVectorFESpace(start, end, cells);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = false;
		final Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		final VectorFESpaceFunction<TPVectorFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		assertTrue("" + ConvergenceOrderEstimator.normL2VecDifference(solut,
		                                                              LaplaceReferenceSolution.vectorReferenceSolution(),
		                                                              grid.generatePlotPoints(20)),
		           ConvergenceOrderEstimator.normL2VecDifference(solut,
		                                                         LaplaceReferenceSolution.vectorReferenceSolution(),
		                                                         grid.generatePlotPoints(20)) < 1e-2);
	}
	
	@Test
	public void testContinuousConvergence()
	{
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int mult = 4;
		final List<Integer> cells = List.of(mult * 2, mult * 3);
		
		final double penalty = 100;
		
		final TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		
		final TPVectorRightHandSideIntegral<ContinuousTPVectorFunction> rightHandSideIntegral =
			new TPVectorRightHandSideIntegral<>(LaplaceReferenceSolution.vectorRightHandSide(),
			                                    TPVectorRightHandSideIntegral.VALUE);
		
		final ArrayList<CellIntegral<TPCell, ContinuousTPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		final ArrayList<FaceIntegral<TPFace, ContinuousTPVectorFunction>> faceIntegrals =
			new ArrayList<>();
		final ArrayList<RightHandSideIntegral<TPCell, ContinuousTPVectorFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		final ArrayList<BoundaryRightHandSideIntegral<TPFace, ContinuousTPVectorFunction>> boundaryFaceIntegrals
			=
			new ArrayList<>();
		
		rightHandSideIntegrals.add(rightHandSideIntegral);
		cellIntegrals.add(gg);
		
		final int polynomialDegree = 2;
		final ContinuousTPFEVectorSpace grid = new ContinuousTPFEVectorSpace(start, end, cells);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setBoundaryValues(LaplaceReferenceSolution.vectorBoundaryValues());
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = false;
		final Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		final VectorFESpaceFunction<ContinuousTPVectorFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(solut,
		                                                         LaplaceReferenceSolution.vectorReferenceSolution(),
		                                                         grid.generatePlotPoints(20)) < 1e-2);
	}
}
