import basic.*;
import distorted.geometry.CircleGrid;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateVector;
import tensorproduct.geometry.CartesianGrid;

public class CircularGrid
{
	public static void main(String[] args)
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = true;
		builder.build();
		CircleGrid circle = new CircleGrid(3,new CoordinateVector(3), 1, 2);
		PlotWindow p = new PlotWindow();
		for(DistortedCell c:circle.cells)
			p.addPlot(new ScalarPlot3D(c.indicatorFunction(), circle.generatePlotPoints(30), 30, "cell"));
		for(DistortedFace c:circle.faces)
			p.addPlot(new ScalarPlot3D(c.indicatorFunction(), circle.generatePlotPoints(100), 101, "face"));
		
	}
}
