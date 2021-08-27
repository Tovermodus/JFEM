package distorted;

import basic.Overlay;
import basic.ScalarPlot2D;
import linalg.CoordinateVector;
import linalg.Matrix;

import java.awt.*;
import java.util.List;

public class DistortedOverlay extends Overlay
{
	private final List<CoordinateVector> backGroundPoints;
	private final DistortedVectorSpace space;
	private final Matrix displacementHistory;
	private final int pointsPerCell;
	CoordinateVector max;
	CoordinateVector min;
	
	public DistortedOverlay(final List<CoordinateVector> backGroundPoints, final DistortedVectorSpace space, final Matrix displacementHistory, final int pointsPerCell)
	{
		this.backGroundPoints = backGroundPoints;
		this.space = space;
		this.displacementHistory = displacementHistory;
		System.out.println(displacementHistory);
		max = ScalarPlot2D.getMaxCoordinates(backGroundPoints);
		min = ScalarPlot2D.getMinCoordinates(backGroundPoints);
		this.pointsPerCell = pointsPerCell * space.getCells().size();
	}
	
	public static String title()
	{
		return "distorted";
	}
	
	@Override
	public void draw(final Graphics g, final int width, final int height, final double slider)
	{
		final int time = (int) (displacementHistory.getRows() * slider * 0.999) % displacementHistory.getRows();
		final DistortedVectorFESpaceFunction X = new DistortedVectorFESpaceFunction(space.getShapeFunctions(),
		                                                                            displacementHistory.getRow(
			                                                                            time));
		int i = 0;
		for (final CoordinateVector c : space.generatePlotPoints(pointsPerCell))
		{
			ScalarPlot2D.drawSinglePoint(g, width, height, X.value(c), (i++ % pointsPerCell*3.3), 0, 100,
			                             min, max, 2, 2);
		}
	}
}
