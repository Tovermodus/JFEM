import basic.*;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class ContinuousVectorLaplace
{
	public static void main(final String[] args)
	{
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int mult = 8;
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
		final ArrayList<BoundaryRightHandSideIntegral<TPFace, ContinuousTPVectorFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		rightHandSideIntegrals.add(rightHandSideIntegral);
		cellIntegrals.add(gg);
		
		final int polynomialDegree = 2;
		final ContinuousTPFEVectorSpace grid = new ContinuousTPFEVectorSpace(start, end, cells);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.setBoundaryValues(LaplaceReferenceSolution.vectorBoundaryValues());
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		for (final Integer s : grid.getFixedNodeIndices())
			System.out.println(grid.getShapeFunctions().get(s).getNodeFunctionalPoint().absMaxElement());
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = false;
		final Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		final VectorFESpaceFunction<ContinuousTPVectorFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		final PlotWindow p = new PlotWindow();
		p.addPlot(new ScalarPlot2D(solut.getComponentFunction(0), grid.generatePlotPoints(30), 30));
		p.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.vectorReferenceSolution().getComponentFunction(0),
		                           grid.generatePlotPoints(30), 30));
		System.out.println(ConvergenceOrderEstimator.normL2VecDifference(solut,
		                                                                 LaplaceReferenceSolution.vectorReferenceSolution(),
		                                                                 grid.generatePlotPoints(20)));
	}
}
