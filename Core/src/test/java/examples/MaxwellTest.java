package examples;

import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.BlockSparseMatrix;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import mixed.*;
import org.junit.Test;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class MaxwellTest
{
	@Test(timeout = 200000)
	public void testQkQkConvergence()
	{
		final CoordinateVector start = CoordinateVector.fromValues(0, 0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1, 1);
		final int polynomialDegree = 3;
		final QkQkSpace grid = new QkQkSpace(start, end,
		                                     Ints.asList(4, 4, 4));
		final TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegralViaReferenceCell<>(TPVectorCellIntegral.ROT_ROT);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vv
			=
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		final List<CellIntegral<TPCell, QkQkFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		final List<FaceIntegral<TPFace, QkQkFunction>> faceIntegrals = new ArrayList<>();
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(MaxwellReferenceSolution.rightHandSide(),
				                                    TPVectorRightHandSideIntegral.VALUE));
		final List<RightHandSideIntegral<TPCell, QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		final List<BoundaryRightHandSideIntegral<TPFace, QkQkFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("velocity Boundary");
		grid.setVelocityBoundaryValues(MaxwellReferenceSolution.vectorBoundaryValues());
		System.out.println("pressure Boundary");
		grid.setPressureBoundaryValues(MaxwellReferenceSolution.pressureBoundaryValues());
		final Stopwatch s = Stopwatch.createStarted();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println(s.elapsed());
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		System.out.println(
			"solve system: " + grid.getSystemMatrix()
			                       .getRows() + "×" + grid.getSystemMatrix()
			                                              .getCols());
		final IterativeSolver i = new IterativeSolver(false);
		i.showProgress = false;
		final Vector solution1 = i.solveGMRES(new BlockSparseMatrix(grid.getSystemMatrix(), 5), grid.getRhs(),
		                                      1e-6);
		System.out.println("solved");
		final MixedTPFESpaceFunction<QkQkFunction> solut =
			new MixedTPFESpaceFunction<>(
				grid.getShapeFunctionMap(), solution1);
		assertTrue(ConvergenceOrderEstimator.normL2Difference(solut.getPressureFunction(),
		                                                      MaxwellReferenceSolution.pressureReferenceSolution(),
		                                                      grid.generatePlotPoints(20)) < 2e-1);
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(solut.getVelocityFunction(),
		                                                         MaxwellReferenceSolution.velocityReferenceSolution(),
		                                                         grid.generatePlotPoints(20)) < 1e-1);
	}
}
