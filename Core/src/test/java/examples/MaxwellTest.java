package examples;

import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.SparseMatrix;
import linalg.Vector;
import mixed.*;
import org.junit.jupiter.api.Test;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class MaxwellTest
{
	@Test
	public void testQkQkConvergence()
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		CoordinateVector start = CoordinateVector.fromValues(0, 0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1, 1);
		int polynomialDegree = 3;
		QkQkSpace grid = new QkQkSpace(start, end,
			Ints.asList(4, 4, 4), polynomialDegree);
		TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegralViaReferenceCell<>(TPVectorCellIntegral.ROT_ROT);
		MixedCellIntegral<TPCell, TPFace, TPEdge, ContinuousTPShapeFunction, ContinuousTPVectorFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, TPFace, TPEdge, ContinuousTPShapeFunction, ContinuousTPVectorFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell, MixedShapeFunction<TPCell, TPFace, TPEdge, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral<TPFace, MixedShapeFunction<TPCell, TPFace, TPEdge, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> faceIntegrals = new ArrayList<>();
		MixedRightHandSideIntegral<TPCell, TPFace, TPEdge, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<ContinuousTPVectorFunction>(MaxwellReferenceSolution.rightHandSide(),
					TPVectorRightHandSideIntegral.VALUE));
		List<RightHandSideIntegral<TPCell, MixedShapeFunction<TPCell, TPFace, TPEdge,
			ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		List<BoundaryRightHandSideIntegral<TPFace, MixedShapeFunction<TPCell, TPFace,
			TPEdge, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> boundaryFaceIntegrals =
			new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("velocity Boundary");
		grid.setVelocityBoundaryValues(MaxwellReferenceSolution.vectorBoundaryValues());
		System.out.println("pressure Boundary");
		grid.setPressureBoundaryValues(MaxwellReferenceSolution.pressureBoundaryValues());
		Stopwatch s = Stopwatch.createStarted();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println(s.elapsed());
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		IterativeSolver i = new IterativeSolver();
		i.showProgress = false;
		Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-6);
		System.out.println("solved");
		MixedFESpaceFunction<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction,ContinuousTPVectorFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		assertTrue(ConvergenceOrderEstimator.normL2Difference(solut.getPressureFunction(),
			MaxwellReferenceSolution.pressureReferenceSolution(), grid.generatePlotPoints(20)) < 2e-1);
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(solut.getVelocityFunction(),
			MaxwellReferenceSolution.velocityReferenceSolution(), grid.generatePlotPoints(20)) < 1e-1);
	}
}
