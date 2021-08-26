package distorted;

import basic.Overlay;
import basic.ScalarPlot2D;
import basic.ScalarPlot3D;
import linalg.CoordinateVector;

import java.awt.*;

public class CircleOverlay extends Overlay
{
	CoordinateVector mins;
	CoordinateVector maxs;
	double maxValue;
	double minValue;
	double radius;
	CoordinateVector center;
	
	public CircleOverlay(final CircleGridSpace<?, ?, ?, ?> circle, final ScalarPlot2D plot)
	{
		radius = circle.radius;
		center = circle.center;
		mins = plot.min;
		maxs = plot.max;
		maxValue = plot.maxValue;
		minValue = plot.minValue;
	}
	
	public CircleOverlay(final CircleGridSpace<?, ?, ?, ?> circle, final ScalarPlot3D plot)
	{
		radius = circle.radius;
		center = circle.center;
		mins = plot.min;
		maxs = plot.max;
		maxValue = plot.maxValue;
		minValue = plot.minValue;
	}
	
	@Override
	public void draw(final Graphics g, final int width, final int height, final double slider)
	{
		for (int i = 0; i < 3600; i++)
		{
			final CoordinateVector circleCoordinate = center
				.add(CoordinateVector.getUnitVector(2, 0).mul(Math.cos(Math.PI * i / 1800)).mul(radius))
				.add(CoordinateVector
					     .getUnitVector(2, 1)
					     .mul(Math.sin(Math.PI * i / 1800))
					     .mul(radius));
			ScalarPlot2D.drawSinglePoint(g, width, height, circleCoordinate,
			                             minValue + (maxValue - minValue) * (i % 100) / 100, minValue,
			                             maxValue, mins, maxs, 5, 5);
		}
	}
}
