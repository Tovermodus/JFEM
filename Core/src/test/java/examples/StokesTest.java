package examples;

import basic.*;
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

public class StokesTest
{
	@Test
	public void testConvergence()
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		CoordinateVector start = CoordinateVector.fromValues(0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 3;
		TaylorHoodSpace grid = new TaylorHoodSpace(start, end,
			Ints.asList(7,5), polynomialDegree);
		TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
				TPVectorCellIntegral.GRAD_GRAD);
		MixedCellIntegral<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction, ContinuousTPVectorFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction, ContinuousTPVectorFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(gradGrad);
		List<CellIntegral<TPCell,  MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		
		List<FaceIntegral< TPFace, MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> faceIntegrals = new ArrayList<>();
		
		MixedRightHandSideIntegral<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<ContinuousTPVectorFunction>(StokesReferenceSolution.rightHandSide(), TPVectorRightHandSideIntegral.VALUE));
		
		List<RightHandSideIntegral<TPCell,  MixedShapeFunction<TPCell, TPFace,TPEdge,
			ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		List<BoundaryRightHandSideIntegral< TPFace, MixedShapeFunction<TPCell, TPFace,TPEdge,
			ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> boundaryFaceIntegrals = new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
		System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		IterativeSolver<SparseMatrix> i = new IterativeSolver<>();
		i.showProgress = false;
		Vector solution1 = i.solveCG(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		System.out.println("solved");
		
		MixedFESpaceFunction<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction,ContinuousTPVectorFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		assertTrue(ConvergenceOrderEstimator.normL20Difference(solut.getPressureFunction(),
			StokesReferenceSolution.pressureReferenceSolution(), grid.generatePlotPoints(20)) < 1e-2);
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(solut.getVelocityFunction(),
			StokesReferenceSolution.velocityReferenceSolution(), grid.generatePlotPoints(20)) < 1e-2);
	}
}
