import basic.PerformanceArguments;
import basic.PlotWindow;
import basic.ScalarPlot2D;
import distorted.geometry.CircleGrid;
import distorted.geometry.DistortedCell;
import linalg.CoordinateVector;

public class CircularGrid
{
	public static void main(final String[] args)
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		builder.build();
		final CircleGrid circle = new CircleGrid(2, new CoordinateVector(2), 2 * Math.pow(2, 1. / 2), 2);
		final PlotWindow p = new PlotWindow();
		System.out.println(circle.faces.size());
		System.out.println(circle.faces);
		for (final DistortedCell c : circle.cells)
			p.addPlot(new ScalarPlot2D(c.indicatorFunction(), circle.generatePlotPoints(120), 150, "cell"));
		final int i = 0;
//		for (final DistortedFace c : circle.faces)
//		{
//			p.addPlot(new ScalarPlot2D(c.indicatorFunction(), circle.generatePlotPoints(30), 50,
//			                           "face " + i++));
//		}
	}
}
