import basic.*;
import distorted.CircleSpace;
import distorted.DistortedCellIntegral;
import distorted.DistortedRightHandSideIntegral;
import distorted.DistortedShapeFunction;
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
		final CircleSpace circle = new CircleSpace(new CoordinateVector(2), 1, 3);
		
		System.out.println("Cells done");
		circle.assembleCells();
		System.out.println("Cells done");
		circle.assembleFunctions(1);
		System.out.println("Cells done");
		circle.initializeSystemMatrix();
		System.out.println("Cells done");
		circle.initializeRhs();
		System.out.println(circle.getShapeFunctionMap()
		                         .size());
		System.out.println("System Initialized");
		circle.evaluateCellIntegrals(List.of(gradGrad), List.of(source));
		System.out.println("System Initialized");
		circle.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		System.out.println("System Initialized");
		circle.setBoundaryValues(ScalarFunction.constantFunction(-1),
		                         face -> face.center()
		                                     .at(0) < 5);
		System.out.println("System Initialized");
		System.out.println(circle.getShapeFunctionMap()
		                         .size());
		
		System.out.println("System Filled");
//		final Vector solution = circle.getSystemMatrix().solve(circle.getRhs());
		final IterativeSolver iterativeSolver = new IterativeSolver();
		iterativeSolver.showProgress = false;
		final Vector solution = iterativeSolver.solveGMRES(circle.getSystemMatrix(), circle.getRhs(), 1e-10);
		
		final ScalarFESpaceFunction<DistortedShapeFunction> solutionFunction =
			new ScalarFESpaceFunction<>(circle.getShapeFunctionMap(), solution);
		
		final int plotploints = 6;
		PlotWindow.addPlot(new ScalarPlot2D(solutionFunction,
		                                    circle.generatePlotPoints(plotploints * circle.getCells()
		                                                                                  .size()),
		                                    5));
		PlotWindow.addPlot(new ScalarPlot2D(solutionFunction,
		                                    circle.generatePlotPoints(plotploints * circle.getCells()
		                                                                                  .size()),
		                                    10));
		PlotWindow.addPlot(new ScalarPlot2D(solutionFunction,
		                                    circle.generatePlotPoints(plotploints * circle.getCells()
		                                                                                  .size()),
		                                    20));
		PlotWindow.addPlot(new ScalarPlot2D(solutionFunction,
		                                    circle.generatePlotPoints(plotploints * circle.getCells()
		                                                                                  .size()),
		                                    30));
		PlotWindow.addPlot(new ScalarPlot2D(solutionFunction,
		                                    circle.generatePlotPoints(plotploints * circle.getCells()
		                                                                                  .size()),
		                                    50));
		PlotWindow.addPlot(new ScalarPlot2D(solutionFunction,
		                                    circle.generatePlotPoints(plotploints * circle.getCells()
		                                                                                  .size()),
		                                    80));
		PlotWindow.addPlot(new ScalarPlot2D(solutionFunction,
		                                    circle.generatePlotPoints(plotploints * circle.getCells()
		                                                                                  .size()),
		                                    110));
		PlotWindow.addPlot(new ScalarPlot2D(solutionFunction,
		                                    circle.generatePlotPoints(plotploints * circle.getCells()
		                                                                                  .size()),
		                                    150));
	}
}
