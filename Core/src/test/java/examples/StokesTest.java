package examples;

import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import mixed.*;
import org.junit.jupiter.api.Test;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class StokesTest
{
	@Test
	public void testConvergence()
	{
		CoordinateVector start = CoordinateVector.fromValues(0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 3;
		TaylorHoodSpace grid = new TaylorHoodSpace(start, end,
			Ints.asList(7,5));
		TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
				TPVectorCellIntegral.GRAD_GRAD);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(gradGrad);
		List<CellIntegral<TPCell,  QkQkFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		
		List<FaceIntegral< TPFace, QkQkFunction>> faceIntegrals = new ArrayList<>();
		
		MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction,QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<ContinuousTPVectorFunction>(StokesReferenceSolution.rightHandSide(), TPVectorRightHandSideIntegral.VALUE));
		
		List<RightHandSideIntegral<TPCell, QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		List<BoundaryRightHandSideIntegral< TPFace, QkQkFunction>> boundaryFaceIntegrals = new ArrayList<>();
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
		IterativeSolver i = new IterativeSolver();
		i.showProgress = false;
		Vector solution1 = i.solveCG(grid.getSystemMatrix(), grid.getRhs(), 1e-10);
		System.out.println("solved");
		
		MixedFESpaceFunction<QkQkFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		System.out.println(ConvergenceOrderEstimator.normL20Difference(solut.getPressureFunction(),
			StokesReferenceSolution.pressureReferenceSolution(), grid.generatePlotPoints(20)) );
		assertTrue(ConvergenceOrderEstimator.normL20Difference(solut.getPressureFunction(),
			StokesReferenceSolution.pressureReferenceSolution(), grid.generatePlotPoints(20)) < 2e-1);
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(solut.getVelocityFunction(),
			StokesReferenceSolution.velocityReferenceSolution(), grid.generatePlotPoints(20)) < 1e-2);
	}
}
