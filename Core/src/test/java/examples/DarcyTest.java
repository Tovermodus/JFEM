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

public class DarcyTest
{
	@Test
	public void testConvergence()
	{
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 2;
		MixedRTSpace grid = new MixedRTSpace(start, end,
			Ints.asList(5,7));
		TPVectorCellIntegral<RTShapeFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		MixedCellIntegral<TPCell,  ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction>
			divValue = new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell,  ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell,  RTMixedFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral< TPFace, RTMixedFunction>> faceIntegrals = new ArrayList<>();
		RightHandSideIntegral<TPCell, RTMixedFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromPressureIntegral(
				new TPRightHandSideIntegral<>(
					LaplaceReferenceSolution.scalarRightHandSide(),TPRightHandSideIntegral.VALUE ));
		List<RightHandSideIntegral<TPCell, RTMixedFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		List<BoundaryRightHandSideIntegral< TPFace, RTMixedFunction>> boundaryFaceIntegrals = new ArrayList<>();
		MixedBoundaryRightHandSideIntegral<TPFace, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction> dirichlet =
			MixedBoundaryRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorBoundaryFaceIntegral<>(LaplaceReferenceSolution.scalarBoundaryValues(),
					TPVectorBoundaryFaceIntegral.NORMAL_VALUE));
		boundaryFaceIntegrals.add(dirichlet);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		IterativeSolver i = new IterativeSolver();
		i.showProgress = false;
		Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-10);
		System.out.println("solved");
		MixedFESpaceFunction<RTMixedFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
//		PlotWindow p = new PlotWindow();
//		p.addPlot(new MixedPlot2D(solut, grid.generatePlotPoints(20),20));
//		p.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.scalarReferenceSolution(), grid.generatePlotPoints(20),20));
		System.out.println(ConvergenceOrderEstimator.normL2Difference(solut.getPressureFunction(),
			LaplaceReferenceSolution.scalarReferenceSolution(),grid.generatePlotPoints(20)));
		try
		{
			Thread.sleep(00000);
		} catch (InterruptedException e)
		{
			e.printStackTrace();
		}
		assertTrue(ConvergenceOrderEstimator.normL2Difference(solut.getPressureFunction(),
			LaplaceReferenceSolution.scalarReferenceSolution(),grid.generatePlotPoints(20))<1e-2);
	}
}
