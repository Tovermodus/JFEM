package basic;

import linalg.CoordinateVector;

import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;

public class ScalarPlot3D implements Plot
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
	public ScalarPlot3D(ScalarFunction function, List<CoordinateVector> points, int pointsPerDimension)
	{
		this(function.valuesInPoints(points), pointsPerDimension, "Unnamed");
	}
	public ScalarPlot3D(ScalarFunction function, List<CoordinateVector> points, int pointsPerDimension,
	                    String title)
	{
		this(function.valuesInPoints(points), pointsPerDimension, title);
	}
	public ScalarPlot3D(Map<CoordinateVector, Double> values, int pointsPerDimension,
	                    String title)
	{
		for (CoordinateVector c:values.keySet())
			if(c.getLength() != 3)
				throw new IllegalArgumentException("Not 3D");
		this.values = values;
		this.pointsPerDimension = pointsPerDimension;
		maxValue = values.values().stream().parallel().max(Double::compare).orElse(0.0);
		minValue = values.values().stream().parallel().min(Double::compare).orElse(0.0);
		max = ScalarPlot2D.getMaxCoordinates(values.keySet());
		min = ScalarPlot2D.getMinCoordinates(values.keySet());
		this.title = title;
	}
	public List<CoordinateVector> getClosestPoints(double slider)
	{
		double currentZ = min.z() + slider*(max.z() - min.z());
		List<CoordinateVector> drawnVectors = new ArrayList<>(values.keySet());
		
		ToDoubleFunction<CoordinateVector> sortFunction = value -> -Math.abs(value.z() - currentZ);
		Comparator<CoordinateVector> comparator = Comparator.comparingDouble(sortFunction)
			.thenComparing(CoordinateVector::x).thenComparing(CoordinateVector::y);
		drawnVectors.sort(comparator);
		double multipleToSecureNoWhitespace = 2;
		drawnVectors =
			drawnVectors.subList(drawnVectors.size()-(int)(multipleToSecureNoWhitespace
					*pointsPerDimension*pointsPerDimension),
				drawnVectors.size()-1);
		return drawnVectors;
	}
	@Override
	public void drawValues(Graphics g, int width, int height, double slider)
	{
		int pixelWidth = (width-150)/pointsPerDimension+1;
		int pixelHeight = (height-150)/pointsPerDimension+1;
		
		for(CoordinateVector co:getClosestPoints(slider))
		{
			ScalarPlot2D.drawSinglePoint(g, width, height, co, values.get(co), minValue, maxValue, min, max,
				pixelWidth,
				pixelHeight);
		}
	}
	
	@Override
	public String title()
	{
		return "3D Plot " + title;
	}
}
