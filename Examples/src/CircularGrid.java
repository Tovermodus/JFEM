import basic.PerformanceArguments;
import basic.PlotWindow;
import basic.ScalarPlot2D;
import distorted.DistortedShapeFunction;
import distorted.geometry.CircleGrid;
import linalg.CoordinateVector;

import java.util.HashSet;

public class CircularGrid
{
	public static void main(final String[] args)
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = true;
		builder.build();
		final CircleGrid circle = new CircleGrid(2, new CoordinateVector(2), 2 * Math.pow(2, 1. / 2), 3);
		final PlotWindow p = new PlotWindow();
		System.out.println(circle.faces.size());
		System.out.println(circle.cells.size());
		//System.out.println(circle.faces);
		final HashSet<DistortedShapeFunction> shapeFunctionSet = new HashSet<>();
		circle.cells.stream().parallel().forEach(c ->
		                                         {
			                                         for (int i = 0; i < 4; i++)
				                                         shapeFunctionSet.add(
					                                         new DistortedShapeFunction(c, 1, i));
		                                         });
		int i = 0;
		for (final DistortedShapeFunction sf : shapeFunctionSet)
			p.addPlot(new ScalarPlot2D(sf, circle.generatePlotPoints(100), 100, "func" + (i++)));
//			p.addPlot(new ScalarPlot2D(c.indicatorFunction(), circle.generatePlotPoints(200), 300, "cell"));
//		int i = 0;
//		for (final DistortedFace c : circle.faces)
//		{
//			p.addPlot(new ScalarPlot2D(c.indicatorFunction(), circle.generatePlotPoints(200), 300,
//			                           "face " + i++));
//		}
	}
}
