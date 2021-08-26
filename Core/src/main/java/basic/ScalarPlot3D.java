package basic;

import linalg.CoordinateVector;

import java.awt.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.ToDoubleFunction;

public class ScalarPlot3D extends Plot
{
	final public Map<CoordinateVector, Double> values;
	final public int pointsPerDimension;
	final public double minValue;
	final public double maxValue;
	final public CoordinateVector min;
	final public CoordinateVector max;
	public final String title;
	public int pixelWidth;
	public int pixelHeight;
	public double z;
	
	public ScalarPlot3D(final ScalarFunction function, final List<CoordinateVector> points, final int pointsPerDimension)
	{
		this(function.valuesInPoints(points), pointsPerDimension, "Unnamed");
	}
	
	public ScalarPlot3D(final ScalarFunction function, final List<CoordinateVector> points, final int pointsPerDimension, final String title)
	{
		this(function.valuesInPoints(points), pointsPerDimension, title);
	}
	
	public ScalarPlot3D(final Map<CoordinateVector, Double> values, final int pointsPerDimension, final String title)
	{
		for (final CoordinateVector c : values.keySet())
			if (c.getLength() != 3) throw new IllegalArgumentException("Not 3D");
		this.values = values;
		this.pointsPerDimension = pointsPerDimension;
		maxValue = values.values().stream().parallel().max(Double::compare).orElse(0.0);
		minValue = values.values().stream().parallel().min(Double::compare).orElse(0.0);
		max = ScalarPlot2D.getMaxCoordinates(values.keySet());
		min = ScalarPlot2D.getMinCoordinates(values.keySet());
		this.title = title;
	}
	
	public List<CoordinateVector> getClosestPoints(final double slider)
	{
		final double currentZ = min.z() + slider * (max.z() - min.z());
		z = currentZ;
		List<CoordinateVector> drawnVectors = new ArrayList<>(values.keySet());
		
		final ToDoubleFunction<CoordinateVector> sortFunction = value -> -Math.abs(value.z() - currentZ);
		final Comparator<CoordinateVector> comparator = Comparator
			.comparingDouble(sortFunction)
			.thenComparing(CoordinateVector::x)
			.thenComparing(CoordinateVector::y);
		drawnVectors.sort(comparator);
		final double multipleToSecureNoWhitespace = 1;
		drawnVectors = drawnVectors.subList(Math.max(0,
		                                             drawnVectors.size() - (int) (multipleToSecureNoWhitespace * pointsPerDimension * pointsPerDimension)),
		                                    drawnVectors.size() - 1);
		return drawnVectors;
	}
	
	@Override
	public void drawValues(final Graphics g, final int width, final int height, final double slider)
	{
		pixelWidth = (width - 150) / pointsPerDimension + 1;
		pixelHeight = (height - 150) / pointsPerDimension + 1;
		
		for (final CoordinateVector co : getClosestPoints(slider))
		{
			ScalarPlot2D.drawSinglePoint(g, width, height, co, values.get(co), minValue, maxValue, min, max,
			                             pixelWidth, pixelHeight);
		}
	}
	
	@Override
	public String title()
	{
		return "3D Plot " + title + " blue is min with " + minValue + " red is max with " + maxValue;
	}
}
