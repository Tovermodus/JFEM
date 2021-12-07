package distorted;

import basic.Overlay;
import basic.ScalarPlot2D;
import distorted.geometry.DistortedCell;
import io.vavr.Tuple2;
import linalg.CoordinateVector;
import linalg.Matrix;

import java.awt.*;
import java.util.List;

public class DistortedOverlay
	extends Overlay
{
	private final CircleVectorSpace space;
	private final Matrix displacementHistory;
	private final int pointsPerCell;
	CoordinateVector max;
	CoordinateVector min;
	int prevTime = -1000;
	List<Tuple2<DistortedCell, CoordinateVector>> refPoints;
	final DistortedVectorFESpaceFunction X;
	
	public DistortedOverlay(final List<CoordinateVector> backGroundPoints,
	                        final CircleVectorSpace space,
	                        final Matrix displacementHistory,
	                        final int pointsPerCell)
	{
		if (displacementHistory.getCols() != space.getShapeFunctions()
		                                          .size())
			throw new IllegalArgumentException("displacement history has wrong dimensions");
		this.space = space;
		this.displacementHistory = displacementHistory;
		//System.out.println("displacemantHist " + displacementHistory);
		max = ScalarPlot2D.getMaxCoordinates(backGroundPoints);
		min = ScalarPlot2D.getMinCoordinates(backGroundPoints);
		this.pointsPerCell = pointsPerCell * space.getCells()
		                                          .size();
		X = new DistortedVectorFESpaceFunction(space.getShapeFunctions(),
		                                       displacementHistory.getRow(0));
		refPoints = space.generateReferencePlotPoints(this.pointsPerCell);
	}
	
	public static String title()
	{
		return "distorted";
	}
	
	@Override
	public void draw(final Graphics g, final int width, final int height, final double slider)
	{
		final int pixelSize = Math.min(width, height) / (int) Math.sqrt(refPoints.size()) / 8;
		final int time =
			(int) (displacementHistory.getRows() * slider * 0.999) % displacementHistory.getRows();
		if (time != prevTime)
		{
			X.resetCoefficients(space.getShapeFunctions(), displacementHistory.getRow(time));
			prevTime = time;
		}
		int i = 0;
		for (final Tuple2<DistortedCell, CoordinateVector> cp : refPoints)
		{
			ScalarPlot2D.drawSinglePointGreen(g, width, height,
			                                  cp._1.transformFromReferenceCell(cp._2)
			                                       .add(X.valueOnReferenceCell(cp._2,
			                                                                   cp._1)),
			                                  (i++ % pointsPerCell * 3.3) % 100, 0,
			                                  100, min, max, pixelSize, pixelSize);
		}
	}
}
