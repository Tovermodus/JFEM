package basic;

import com.google.common.collect.Iterables;
import io.vavr.Tuple2;
import linalg.CoordinateVector;

import java.awt.*;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Optional;

public class ScalarPlot2D
	extends Plot
{
	final public Map<CoordinateVector, Double> values;
	final public int pointsPerDimension;
	final public double minValue;
	final public double maxValue;
	final public CoordinateVector min;
	final public CoordinateVector max;
	public int pixelWidth;
	public int pixelHeight;
	public final String title;
	
	public ScalarPlot2D(final ScalarFunction function,
	                    final List<CoordinateVector> points,
	                    final int pointsPerDimension)
	{
		this(function.valuesInPoints(points), pointsPerDimension, "Unnamed");
	}
	
	public ScalarPlot2D(final ScalarFunction function,
	                    final List<CoordinateVector> points,
	                    final int pointsPerDimension,
	                    final String title)
	{
		this(function.valuesInPoints(points), pointsPerDimension, title);
	}
	
	public ScalarPlot2D(final Map<CoordinateVector, Double> values, final int pointsPerDimension,
	                    final String title)
	{
		for (final CoordinateVector c : values.keySet())
			if (c.getLength() != 2)
				throw new IllegalArgumentException("Not 2D");
		this.values = values;
		this.pointsPerDimension = pointsPerDimension;
		maxValue = values.values()
		                 .stream()
		                 .parallel()
		                 .max(Double::compare)
		                 .orElse(0.0);
		minValue = values.values()
		                 .stream()
		                 .parallel()
		                 .min(Double::compare)
		                 .orElse(0.0);
		max = getMaxCoordinates(values.keySet());
		min = getMinCoordinates(values.keySet());
		this.title = title;
	}
	
	@Override
	public void drawValues(final Graphics g, final int width, final int height, final double slider)
	{
		pixelWidth = (width - 150) / pointsPerDimension + 1;
		pixelHeight = (height - 150) / pointsPerDimension + 1;
		for (final Map.Entry<CoordinateVector, Double> entry : values.entrySet())
		{
			drawSinglePoint(g, width, height, entry.getKey(), entry.getValue(), minValue, maxValue, min,
			                max,
			                pixelWidth,
			                pixelHeight);
		}
		g.setColor(Color.black);
		g.drawString("" + min.x(), 50, height - 35);
		g.drawLine(75, height - 40, width - 100, height - 40);
		g.drawString("" + max.x(), width - 75, height - 35);
		g.drawString("" + min.y(), 15, height - 55);
		g.drawLine(25, 85, 25, height - 75);
		g.drawString("" + max.y(), 15, 65);
	}
	
	public static CoordinateVector getMinCoordinates(final Collection<CoordinateVector> coordinates)
	{
		final int dim = Iterables.getLast(coordinates)
		                         .getLength();
		final CoordinateVector ret = new CoordinateVector(dim);
		for (int i = 0; i < dim; i++)
		{
			final int finalI = i;
			final Optional<Double> min =
				coordinates.stream()
				           .parallel()
				           .map(c -> c.at(finalI))
				           .min(Double::compare);
			ret.set(min.orElse(0.0), i);
		}
		return ret;
	}
	
	public static CoordinateVector getMaxCoordinates(final Collection<CoordinateVector> coordinates)
	{
		final int dim = Iterables.getLast(coordinates)
		                         .getLength();
		final CoordinateVector ret = new CoordinateVector(dim);
		for (int i = 0; i < dim; i++)
		{
			final int finalI = i;
			final Optional<Double> min =
				coordinates.stream()
				           .parallel()
				           .map(c -> c.at(finalI))
				           .max(Double::compare);
			ret.set(min.orElse(0.0), i);
		}
		return ret;
	}
	
	public static void drawSinglePoint(final Graphics g,
	                                   final int width,
	                                   final int height,
	                                   final CoordinateVector coord,
	                                   final double val,
	                                   final double minValue,
	                                   final double maxValue,
	                                   final CoordinateVector mins,
	                                   final CoordinateVector maxs,
	                                   final int pixelWidth,
	                                   final int pixelHeight)
	{
		double maxValueDifference = maxValue - minValue;
		if (maxValueDifference < 1e-15)
			maxValueDifference = 1;
		final int red = (int) ((val - minValue) / maxValueDifference * 255);
		final int green = 0;
		final int blue = 255 - red;
		final Color c = new Color(red, green, blue);
		g.setColor(c);
		final Tuple2<Integer, Integer> xy = getXY(coord, mins, maxs, width, height, pixelWidth, pixelHeight);
		g.fillRect(xy._1, xy._2, pixelWidth + 2, pixelHeight + 2);
	}
	
	public static void drawSinglePointGreen(final Graphics g,
	                                        final int width,
	                                        final int height,
	                                        final CoordinateVector coord,
	                                        final double val,
	                                        final double minValue,
	                                        final double maxValue,
	                                        final CoordinateVector mins,
	                                        final CoordinateVector maxs,
	                                        final int pixelWidth,
	                                        final int pixelHeight)
	{
		double maxValueDifference = maxValue - minValue;
		if (maxValueDifference < 1e-15)
			maxValueDifference = 1;
		final int green = (int) ((val - minValue) / maxValueDifference * 255);
		final int blue = 0;
		final int red = 0;
		final Color c = new Color(red, green, blue);
		g.setColor(c);
		final Tuple2<Integer, Integer> xy = getXY(coord, mins, maxs, width, height, pixelWidth, pixelHeight);
		g.fillRect(xy._1, xy._2, pixelWidth + 2, pixelHeight + 2);
	}
	
	public static int getX(final double x, final double minX, final double maxX, final int width,
	                       final int drawnObjectWidth)
	{
		return (int) ((x - minX) / (maxX - minX) * (width - 150) + 75 - drawnObjectWidth / 2);
	}
	
	public static int getY(final double y, final double minY, final double maxY, final int height,
	                       final int drawnObjectHeight)
	{
		return (int) (height - ((y - minY) / (maxY - minY) * (height - 150) + 75 +
			                        drawnObjectHeight / 2));
	}
	
	public static Tuple2<Integer, Integer> getXY(final CoordinateVector coord, final CoordinateVector min,
	                                             final CoordinateVector max, final int width, final int height,
	                                             final int drawnObjectWidth, final int drawnObjectHeight)
	{
		return new Tuple2<>(getX(coord.x(), min.x(), max.x(), width, drawnObjectWidth),
		                    getY(coord.y(), min.y(), max.y(), height, drawnObjectHeight));
	}
	
	@Override
	public String title()
	{
		return "2D Plot " + title + " blue is min with " + minValue + " red is max with " + maxValue;
	}
}
