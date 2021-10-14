package distorted;

import basic.Overlay;
import basic.ScalarPlot2D;
import com.google.common.base.Stopwatch;
import linalg.CoordinateVector;
import linalg.Matrix;

import java.awt.*;
import java.util.List;
import java.util.concurrent.TimeUnit;

public class DistortedOverlay extends Overlay
{
	private final List<CoordinateVector> backGroundPoints;
	private final DistortedVectorSpace space;
	private final Matrix displacementHistory;
	private final int pointsPerCell;
	CoordinateVector max;
	CoordinateVector min;
	int prevTime = -1000;
	Iterable<CoordinateVector> points;
	final DistortedVectorFESpaceFunction X;
	
	public DistortedOverlay(final List<CoordinateVector> backGroundPoints, final DistortedVectorSpace space, final Matrix displacementHistory, final int pointsPerCell)
	{
		this.backGroundPoints = backGroundPoints;
		this.space = space;
		this.displacementHistory = displacementHistory;
		System.out.println(displacementHistory);
		max = ScalarPlot2D.getMaxCoordinates(backGroundPoints);
		min = ScalarPlot2D.getMinCoordinates(backGroundPoints);
		this.pointsPerCell = pointsPerCell * space.getCells().size();
		points = space.generatePlotPoints(this.pointsPerCell);
		X = new DistortedVectorFESpaceFunction(space.getShapeFunctions(),
		                                       displacementHistory.getRow(0));
	}
	
	public static String title()
	{
		return "distorted";
	}
	
	@Override
	public void draw(final Graphics g, final int width, final int height, final double slider)
	{
		Stopwatch s = Stopwatch.createStarted();
		final int time =
			(int) (displacementHistory.getRows() * slider * 0.999) % displacementHistory.getRows();
		if (time != prevTime)
			X.resetCoefficients(space.getShapeFunctions(), displacementHistory.getRow(time));
		int i = 0;
		System.out.println("reset" + s.elapsed(TimeUnit.MICROSECONDS));
		s = Stopwatch.createStarted();
		for (final CoordinateVector c : points)
		{
			ScalarPlot2D.drawSinglePoint(g, width, height, X.value(c), (i++ % pointsPerCell * 3.3) % 100, 0,
			                             100, min, max, 2, 2);
		}
		System.out.println("draw" + s.elapsed(TimeUnit.MICROSECONDS));
	}
}
