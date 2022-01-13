package tensorproduct;

import basic.Overlay;
import basic.ScalarPlot2D;
import basic.VectorFESpaceFunction;
import basic.VectorShapeFunction;
import linalg.CoordinateVector;
import linalg.Matrix;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.awt.*;
import java.util.List;

public class CartesianOverlay<ST extends VectorShapeFunction<TPCell, TPFace>>
	extends Overlay
{
	private final CartesianGridSpace<ST, ?, ?, ?> space;
	private final Matrix displacementHistory;
	private final int pointsPerCell;
	CoordinateVector max;
	CoordinateVector min;
	
	public CartesianOverlay(final List<CoordinateVector> backGroundPoints,
	                        final CartesianGridSpace<ST, ?, ?, ?> space, final Matrix displacementHistory,
	                        final int pointsPerCell)
	{
		this.space = space;
		this.displacementHistory = displacementHistory;
		max = ScalarPlot2D.getMaxCoordinates(backGroundPoints);
		min = ScalarPlot2D.getMinCoordinates(backGroundPoints);
		this.pointsPerCell = pointsPerCell * space.getCells()
		                                          .size();
	}
	
	public static String title()
	{
		return "distorted";
	}
	
	@Override
	public void draw(final Graphics g, final int width, final int height, final double slider)
	{
		final int time = (int) (displacementHistory.getRows() * slider * 0.999) % displacementHistory.getRows();
		final VectorFESpaceFunction<ST> X = new VectorFESpaceFunction<>(space.getShapeFunctions(),
		                                                                displacementHistory.getRow(
			                                                                time));
		int i = 0;
		for (final CoordinateVector c : space.generatePlotPoints(pointsPerCell))
		{
			ScalarPlot2D.drawSinglePoint(g, width, height, X.value(c), (i++ % pointsPerCell * 17) % 100, 0,
			                             100, min, max, 2, 2);
		}
	}
}
