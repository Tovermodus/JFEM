package examples;

import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.*;
import org.junit.Test;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class StokesTest
{
	@Test(timeout = 20000)
	public void testBoundaryBeforeAfter()
	{
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 23);
		final int polynomialDegree = 2;
		final TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
			                           TPVectorCellIntegral.GRAD_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vv
			=
			MixedCellIntegral.fromVelocityIntegral(gradGrad);
		
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(StokesReferenceSolution.rightHandSide(),
				                                    TPVectorRightHandSideIntegral.VALUE));
		
		final TaylorHoodSpace grid = new TaylorHoodSpace(start, end,
		                                                 Ints.asList(3, 3));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
		grid.evaluateCellIntegrals(List.of(vv, divValue), List.of(rightHandSideIntegral));
		grid.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		final SparseMatrix mat = grid.getSystemMatrix();
		final DenseVector rhs = grid.getRhs();
		
		final TaylorHoodSpace grid2 = new TaylorHoodSpace(start, end, Ints.asList(3, 3));
		grid2.assembleCells();
		grid2.assembleFunctions(polynomialDegree);
		grid2.initializeSystemMatrix();
		grid2.initializeRhs();
		grid2.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
		System.out.println("Cell Integrals");
		grid2.evaluateCellIntegrals(List.of(vv, divValue), List.of(rightHandSideIntegral));
		System.out.println("Face Integrals");
		grid2.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		final SparseMatrix mat2 = grid2.getSystemMatrix();
		final DenseVector rhs2 = grid2.getRhs();
		assertTrue(mat2.sub(mat).absMaxElement() < 1e-10);
		//assertTrue(mat2.sub(mat2.transpose()).absMaxElement() < 1e-10);
		assertTrue(rhs2.sub(rhs).absMaxElement() < 1e-10);
	}
	
	@Test
	public void testBoundaryConvergence()
	{
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 3;
		final TaylorHoodSpace grid = new TaylorHoodSpace(start, end,
		                                                 Ints.asList(3, 3));
		final TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
			                           TPVectorCellIntegral.GRAD_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vv
			=
			MixedCellIntegral.fromVelocityIntegral(gradGrad);
		final List<CellIntegral<TPCell, QkQkFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		
		final List<FaceIntegral<TPFace, QkQkFunction>> faceIntegrals = new ArrayList<>();
		
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(StokesReferenceSolution.rightHandSide(),
				                                    TPVectorRightHandSideIntegral.VALUE));
		
		final List<RightHandSideIntegral<TPCell, QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		final List<BoundaryRightHandSideIntegral<TPFace, QkQkFunction>> boundaryFaceIntegrals
			= new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		System.out.println(
			"solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		final IterativeSolver i = new IterativeSolver();
		i.showProgress = false;
		final Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-10);
		System.out.println("solved");
		System.out.println(grid.getSystemMatrix().sub(grid.getSystemMatrix().transpose()).absMaxElement());
		final MixedFESpaceFunction<QkQkFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		System.out.println(ConvergenceOrderEstimator.normL20Difference(solut.getPressureFunction(),
		                                                               StokesReferenceSolution.pressureReferenceSolution(),
		                                                               grid.generatePlotPoints(20)));
		assertTrue(ConvergenceOrderEstimator.normL20Difference(solut.getPressureFunction(),
		                                                       StokesReferenceSolution.pressureReferenceSolution(),
		                                                       grid.generatePlotPoints(20)) < 2e-1);
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(solut.getVelocityFunction(),
		                                                         StokesReferenceSolution.velocityReferenceSolution(),
		                                                         grid.generatePlotPoints(20)) < 1e-2);
	}
	
	@Test
	public void testConvergenceBoundary()
	{
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 3;
		final TaylorHoodSpace grid = new TaylorHoodSpace(start, end,
		                                                 Ints.asList(3, 3));
		final TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
			                           TPVectorCellIntegral.GRAD_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vv
			=
			MixedCellIntegral.fromVelocityIntegral(gradGrad);
		final List<CellIntegral<TPCell, QkQkFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		
		final List<FaceIntegral<TPFace, QkQkFunction>> faceIntegrals = new ArrayList<>();
		
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(StokesReferenceSolution.rightHandSide(),
				                                    TPVectorRightHandSideIntegral.VALUE));
		
		final List<RightHandSideIntegral<TPCell, QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		final List<BoundaryRightHandSideIntegral<TPFace, QkQkFunction>> boundaryFaceIntegrals
			= new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		System.out.println(
			"solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		grid.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
		final IterativeSolver i = new IterativeSolver();
		i.showProgress = true;
		final Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-10);
		System.out.println("solved");
		System.out.println(grid.getSystemMatrix().sub(grid.getSystemMatrix().transpose()).absMaxElement());
		
		final MixedFESpaceFunction<QkQkFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		System.out.println(ConvergenceOrderEstimator.normL20Difference(solut.getPressureFunction(),
		                                                               StokesReferenceSolution.pressureReferenceSolution(),
		                                                               grid.generatePlotPoints(20)));
		assertTrue(ConvergenceOrderEstimator.normL20Difference(solut.getPressureFunction(),
		                                                       StokesReferenceSolution.pressureReferenceSolution(),
		                                                       grid.generatePlotPoints(20)) < 2e-1);
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(solut.getVelocityFunction(),
		                                                         StokesReferenceSolution.velocityReferenceSolution(),
		                                                         grid.generatePlotPoints(20)) < 1e-2);
	}
}
