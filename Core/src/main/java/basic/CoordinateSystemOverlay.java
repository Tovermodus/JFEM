package basic;

import io.vavr.Tuple2;
import io.vavr.Tuple4;
import linalg.CoordinateVector;

import java.awt.*;
import java.util.concurrent.CopyOnWriteArrayList;

public class CoordinateSystemOverlay
	extends Overlay
{
	private final CoordinateVector startCoordinates;
	private final CoordinateVector endCoordinates;
	private final CopyOnWriteArrayList<Tuple4<CoordinateVector, String, Integer, Color>> points;
	
	public CoordinateSystemOverlay(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates)
	{
		this.startCoordinates = startCoordinates;
		this.endCoordinates = endCoordinates;
		points = new CopyOnWriteArrayList<>();
	}
	
	public void addPoint(final CoordinateVector newPoint, final String title, final int size, final Color color)
	{
		points.add(new Tuple4<>(newPoint, title, size, color));
	}
	
	@Override
	public void draw(final Graphics g, final int width, final int height, final double slider)
	{
		g.setColor(Color.BLACK);
		final Tuple2<Integer, Integer> min = ScalarPlot2D.getXY(startCoordinates,
		                                                        startCoordinates,
		                                                        endCoordinates,
		                                                        width,
		                                                        height,
		                                                        1,
		                                                        1);
		final Tuple2<Integer, Integer> max = ScalarPlot2D.getXY(endCoordinates,
		                                                        startCoordinates,
		                                                        endCoordinates,
		                                                        width,
		                                                        height,
		                                                        1,
		                                                        1);
		g.drawLine(min._1, min._2, max._1, min._2);
		g.drawLine(min._1, min._2, min._1, max._2);
		for (final Tuple4<CoordinateVector, String, Integer, Color> point : points)
		{
			drawPoint(g, width, height, point._1, point._2, true, point._3, point._4);
		}
//		drawPoint(g, width, height, startCoordinates, String.format("%4.1e", startCoordinates.x()), false);
//		drawPoint(g, width, height, startCoordinates, String.format("%4.1e", startCoordinates.y()), true);
//		drawPoint(g, width, height, endCoordinates, String.format("%4.1e", endCoordinates.x()), false);
//		drawPoint(g, width, height, endCoordinates, String.format("%4.1e", endCoordinates.y()), true);
//		g.drawString(String.format("%4.1e", endCoordinates.x()), max._1 - 12, min._2 + 19);
//		g.drawString(String.format("%4.1e", endCoordinates.y()), min._1 - 52, max._2 - 1);
	}
	
	private void drawPoint(final Graphics g, final int width, final int height, final CoordinateVector point,
	                       final String title,
	                       final boolean labelLeft, final int size, final Color c)
	{
		
		final Tuple2<Integer, Integer> pointOnGrid = ScalarPlot2D.getXY(point,
		                                                                startCoordinates,
		                                                                endCoordinates,
		                                                                width,
		                                                                height,
		                                                                4,
		                                                                4);
		g.setColor(c);
		g.fillOval(pointOnGrid._1 - size / 2, pointOnGrid._2 - size / 2, size, size);
		g.setColor(Color.BLACK);
		if (labelLeft)
		{
			g.drawString(title, pointOnGrid._1 - 7 * title.length(), pointOnGrid._2 + 5);
		} else
		{
			g.drawString(title, pointOnGrid._1 - 7 * title.length() / 2,
			             pointOnGrid._2 + 20);
		}
	}
}
