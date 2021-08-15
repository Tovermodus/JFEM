import basic.PerformanceArguments;
import basic.PlotWindow;
import basic.ScalarFESpaceFunction;
import basic.ScalarPlot2D;
import distorted.DistortedCellIntegral;
import distorted.DistortedRightHandSideIntegral;
import distorted.DistortedShapeFunction;
import distorted.DistortedSpace;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;

import java.util.ArrayList;
import java.util.List;

public class DiskLaplace
{
	public static void main(final String[] args)
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = true;
		builder.build();
		final DistortedCellIntegral gradGrad = new DistortedCellIntegral(1, DistortedCellIntegral.GRAD_GRAD);
		final DistortedRightHandSideIntegral source = new DistortedRightHandSideIntegral(
			LaplaceReferenceSolution.scalarRightHandSide(),
			DistortedRightHandSideIntegral.VALUE);
		final PlotWindow p = new PlotWindow();
		
		for (int i = 0; i < 4; i++)
		{
			final DistortedSpace circle = new DistortedSpace(new CoordinateVector(2), 1, i);
			
			System.out.println("Cells done");
			circle.assembleCells();
			circle.assembleFunctions(2);
			circle.initializeSystemMatrix();
			circle.initializeRhs();
			System.out.println("System Initialized");
			circle.evaluateCellIntegrals(List.of(gradGrad), List.of(source));
			circle.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
			circle.setBoundaryValues(LaplaceReferenceSolution.scalarReferenceSolution());
			System.out.println(circle.getShapeFunctions().size());
			
			System.out.println("System Filled");
//		final Vector solution = circle.getSystemMatrix().solve(circle.getRhs());
			final IterativeSolver iterativeSolver = new IterativeSolver();
			iterativeSolver.showProgress = false;
			final Vector solution = iterativeSolver.solveGMRES(circle.getSystemMatrix(), circle.getRhs(),
			                                                   1e-10);
			
			final ScalarFESpaceFunction<DistortedShapeFunction> solutionFunction =
				new ScalarFESpaceFunction<>(circle.getShapeFunctions(), solution);
			p.addPlot(new ScalarPlot2D(solutionFunction,
			                           circle.generatePlotPoints(6 * circle.getCells().size()),
			                           30));
			p.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.scalarReferenceSolution(),
			                           circle.generatePlotPoints(6 * circle.getCells().size()), 30));
			System.out.println(ConvergenceOrderEstimator.normL2Difference(solutionFunction,
			                                                              LaplaceReferenceSolution.scalarReferenceSolution(),
			                                                              circle.generatePlotPoints(6 *
				                                                                                        circle
					                                                                                        .getCells()
					                                                                                        .size())));
		}
	}
}
