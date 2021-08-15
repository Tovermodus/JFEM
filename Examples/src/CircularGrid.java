import basic.*;
import distorted.DistortedCellIntegral;
import distorted.DistortedRightHandSideIntegral;
import distorted.DistortedShapeFunction;
import distorted.DistortedSpace;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;

import java.util.ArrayList;
import java.util.List;

public class CircularGrid
{
	public static void main(final String[] args)
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = true;
		builder.build();
		final DistortedCellIntegral gradGrad = new DistortedCellIntegral(1, DistortedCellIntegral.GRAD_GRAD);
		final DistortedRightHandSideIntegral source = new DistortedRightHandSideIntegral(
			ScalarFunction.constantFunction(1),
			DistortedRightHandSideIntegral.VALUE);
		
		final DistortedSpace circle = new DistortedSpace(new CoordinateVector(2), 1, 3);
		
		System.out.println("Cells done");
		circle.assembleCells();
		circle.assembleFunctions(3);
		circle.initializeSystemMatrix();
		circle.initializeRhs();
		System.out.println("System Initialized");
		circle.evaluateCellIntegrals(List.of(gradGrad), List.of(source));
		circle.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		circle.setBoundaryValues(ScalarFunction.constantFunction(-1), face -> face.center().at(0) < 5);
		System.out.println(circle.getShapeFunctions().size());
		
		System.out.println("System Filled");
//		final Vector solution = circle.getSystemMatrix().solve(circle.getRhs());
		final IterativeSolver iterativeSolver = new IterativeSolver();
		iterativeSolver.showProgress = false;
		final Vector solution = iterativeSolver.solveGMRES(circle.getSystemMatrix(), circle.getRhs(), 1e-10);
		
		final ScalarFESpaceFunction<DistortedShapeFunction> solutionFunction =
			new ScalarFESpaceFunction<>(circle.getShapeFunctions(), solution);
		
		final PlotWindow p = new PlotWindow();
//		p.addPlot(new MatrixPlot(circle.getSystemMatrix()));
		p.addPlot(new ScalarPlot2D(solutionFunction, circle.generatePlotPoints(6 * circle.getCells().size()),
		                           5));
		p.addPlot(new ScalarPlot2D(solutionFunction, circle.generatePlotPoints(6 * circle.getCells().size()),
		                           20));
		p.addPlot(new ScalarPlot2D(solutionFunction, circle.generatePlotPoints(6 * circle.getCells().size()),
		                           30));
		p.addPlot(new ScalarPlot2D(solutionFunction, circle.generatePlotPoints(6 * circle.getCells().size()),
		                           50));
		p.addPlot(new ScalarPlot2D(solutionFunction, circle.generatePlotPoints(6 * circle.getCells().size()),
		                           80));
		p.addPlot(new ScalarPlot2D(solutionFunction, circle.generatePlotPoints(6 * circle.getCells().size()),
		                           110));
		p.addPlot(new ScalarPlot2D(solutionFunction, circle.generatePlotPoints(6 * circle.getCells().size()),
		                           150));
	}
}
