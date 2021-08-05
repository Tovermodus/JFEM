import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.MatrixPlot;
import linalg.Vector;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class ContinuousVectorLaplace
{
	public static void main(String[] args)
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
		ContinuousTPFEVectorSpace grid = new ContinuousTPFEVectorSpace(start, end, cells);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.setBoundaryValues(LaplaceReferenceSolution.vectorBoundaryValues());
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		for(Integer s: grid.getFixedNodeIndices())
			System.out.println(grid.getShapeFunctions().get(s).getNodeFunctionalPoint().absMaxElement());
		IterativeSolver it = new IterativeSolver();
		it.showProgress = false;
		Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		VectorFESpaceFunction<ContinuousTPVectorFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		PlotWindow p = new PlotWindow();
		p.addPlot(new ScalarPlot2D(solut.getComponentFunction(0),grid.generatePlotPoints(30),30));
		p.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.vectorReferenceSolution().getComponentFunction(0),
			grid.generatePlotPoints(30),30));
		System.out.println(ConvergenceOrderEstimator.normL2VecDifference(solut,
			LaplaceReferenceSolution.vectorReferenceSolution(), grid.generatePlotPoints(20)));
	}
}
