import basic.PerformanceArguments;
import basic.PlotWindow;
import basic.ScalarPlot3D;
import distorted.geometry.CircleGrid;
import distorted.geometry.DistortedFace;
import linalg.CoordinateVector;

public class CircularGrid
{
	public static void main(final String[] args)
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		builder.build();
		final CircleGrid circle = new CircleGrid(3, new CoordinateVector(3), 2 * Math.pow(2, 1. / 3), 2);
		final PlotWindow p = new PlotWindow();
		System.out.println(circle.faces.size());
		//for(DistortedCell c:circle.cells)
		//	p.addPlot(new ScalarPlot3D(c.indicatorFunction(), circle.generatePlotPoints(30), 30, "cell"));
		int i = 0;
		for (final DistortedFace c : circle.faces)
		{
			p.addPlot(new ScalarPlot3D(c.indicatorFunction(), circle.generatePlotPoints(30), 30,
			                           "face " + i++));
		}
	}
}
