package basic;

import com.google.common.collect.Iterables;
import linalg.CoordinateVector;

import java.awt.*;
import java.util.*;
import java.util.List;

public class ScalarPlot2D implements Plot
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
		maxValue = values.values().stream().parallel().max(Double::compare).orElse(0.0);
		minValue = values.values().stream().parallel().min(Double::compare).orElse(0.0);
		max = getMaxCoordinates(values.keySet());
		min = getMinCoordinates(values.keySet());
		this.title = title;
	}
	@Override
	public void drawValues(Graphics g, int width, int height, double slider)
	{
		pixelWidth = (width-150)/pointsPerDimension+1;
		pixelHeight = (height-150)/pointsPerDimension+1;
		for(Map.Entry<CoordinateVector, Double> entry:values.entrySet())
		{
			drawSinglePoint(g, width, height, entry.getKey(), entry.getValue(), minValue, maxValue, min,
				max,
				pixelWidth,
				pixelHeight);
		}
	}
	public static CoordinateVector getMinCoordinates(Collection<CoordinateVector> coordinates)
	{
		int dim = Iterables.getLast(coordinates).getLength();
		CoordinateVector ret = new CoordinateVector(dim);
		for (int i = 0; i < dim; i++)
		{
			int finalI = i;
			Optional<Double> min =
				coordinates.stream().parallel().map(c->c.at(finalI)).min(Double::compare);
			ret.set(min.orElse(0.0),i);
		}
		return ret;
	}
	public static CoordinateVector getMaxCoordinates(Collection<CoordinateVector> coordinates)
	{
		int dim = Iterables.getLast(coordinates).getLength();
		CoordinateVector ret = new CoordinateVector(dim);
		for (int i = 0; i < dim; i++)
		{
			int finalI = i;
			Optional<Double> min =
				coordinates.stream().parallel().map(c->c.at(finalI)).max(Double::compare);
			ret.set(min.orElse(0.0),i);
		}
		return ret;
	}
	public static void drawSinglePoint(Graphics g, int width, int height,
	                                    CoordinateVector coord, double val, double minValue,
	                                    double maxValue, CoordinateVector mins, CoordinateVector maxs,
	                                    int pixelWidth, int pixelHeight)
	{
		
		double maxValueDifference = maxValue - minValue;
		int red = (int)((val - minValue)/maxValueDifference*255);
		int green = 0;
		int blue = (int)((maxValue - val)/maxValueDifference*255);
		Color c = new Color(red,green,blue);
		int x = (int)((coord.x() - mins.x())/(maxs.x()-mins.x())*(width -150)+75 - pixelWidth/2);
		int y =
			(int)(height - ((coord.y() - mins.y())/(maxs.y()-mins.y())*(height -150)+75 + pixelHeight/2));
		g.setColor(c);
		g.fillRect(x,y,pixelWidth,pixelHeight);
	}
	
	@Override
	public String title()
	{
		return "2D Plot " + title;
	}
}
