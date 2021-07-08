package basic;

import linalg.CoordinateVector;

import java.awt.*;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;

public class ScalarPlot2D implements Plot
{
	final public Map<CoordinateVector, Double> values;
	final public int pointsPerDimension;
	final public double minValue;
	final public double maxValue;
	final public double minX;
	final public double minY;
	final public double maxX;
	final public double maxY;
	public int pixelWidth;
	public int pixelHeight;
	public double maxValueDifference;
	public final String title;
	public ScalarPlot2D(ScalarFunction function, List<CoordinateVector> points, int pointsPerDimension)
	{
		this(function.valuesInPoints(points), pointsPerDimension, "Unnamed");
	}
	public ScalarPlot2D(ScalarFunction function, List<CoordinateVector> points, int pointsPerDimension,
	                    String title)
	{
		this(function.valuesInPoints(points), pointsPerDimension, title);
	}
	public ScalarPlot2D(Map<CoordinateVector, Double> values, int pointsPerDimension,
	                    String title)
	{
		for (CoordinateVector c:values.keySet())
			if(c.getLength() != 2)
				throw new IllegalArgumentException("Not 2D");
		this.values = values;
		this.pointsPerDimension = pointsPerDimension;
		Optional<Double> maV = values.values().stream().parallel().max(Double::compare);
		Optional<Double> miV = values.values().stream().parallel().min(Double::compare);
		Optional<CoordinateVector> maX =
			values.keySet().stream().parallel().max(Comparator.comparingDouble(CoordinateVector::x));
		Optional<CoordinateVector> miX =
			values.keySet().stream().parallel().min(Comparator.comparingDouble(CoordinateVector::x));
		Optional<CoordinateVector> maY =
			values.keySet().stream().parallel().max(Comparator.comparingDouble(CoordinateVector::y));
		Optional<CoordinateVector> miY =
			values.keySet().stream().parallel().min(Comparator.comparingDouble(CoordinateVector::y));
		maxValue = maV.orElse(0.0);
		minValue = miV.orElse(0.0);
		maxX = maX.orElse(new CoordinateVector(2)).x();
		minX = miX.orElse(new CoordinateVector(2)).x();
		maxY = maY.orElse(new CoordinateVector(2)).y();
		minY = miY.orElse(new CoordinateVector(2)).y();
		this.title = title;
	}
	
	@Override
	public void drawValues(Graphics g, int width, int height, double slider)
	{
		pixelWidth = (width-150)/pointsPerDimension+1;
		pixelHeight = (height-150)/pointsPerDimension+1;
		maxValueDifference = maxValue - minValue;
		for(Map.Entry<CoordinateVector, Double> entry:values.entrySet())
		{
			int red = (int)((entry.getValue() - minValue)/maxValueDifference*255);
			int green = 0;
			int blue = (int)((maxValue - entry.getValue())/maxValueDifference*255);
			Color c = new Color(red,green,blue);
			int x = (int)((entry.getKey().x() - minX)/(maxX-minX)*(width-150)+75 - pixelWidth/2);
			int y =
				(int)(height - ((entry.getKey().y() - minY)/(maxY-minY)*(height-150)+75 + pixelHeight/2));
			g.setColor(c);
			g.fillRect(x,y,pixelWidth,pixelHeight);
		}
	}
	
	@Override
	public String title()
	{
		return "2D Plot " + title;
	}
}
