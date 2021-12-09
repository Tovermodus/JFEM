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
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = true;
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int mult = 4;
		final List<Integer> cells = List.of(mult * 2, mult * 2);
		
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
		
		final int polynomialDegree = 1;
		final ContinuousTPFEVectorSpace grid = new ContinuousTPFEVectorSpace(start, end, cells);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setBoundaryValues(LaplaceReferenceSolution.vectorBoundaryValues());
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = true;
		final Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
		final VectorFESpaceFunction<ContinuousTPVectorFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		PlotWindow.addPlot(new ScalarPlot2D(solut.getComponentFunction(0), grid.generatePlotPoints(30), 30));
		PlotWindow.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.vectorReferenceSolution()
		                                                            .getComponentFunction(0),
		                                    grid.generatePlotPoints(30), 30));
		System.out.println(ConvergenceOrderEstimator.normL2VecDifference(solut,
		                                                                 LaplaceReferenceSolution.vectorReferenceSolution(),
		                                                                 grid.generatePlotPoints(20)));
	}
}
