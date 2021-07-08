package basic;

import linalg.CoordinateVector;

import java.awt.*;
import java.util.*;
import java.util.List;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;

public class ScalarPlot3D implements Plot
{
	final Map<CoordinateVector, Double> values;
	final int pointsPerDimension;
	final double minValue;
	final double maxValue;
	final double minX;
	final double minY;
	final double minZ;
	final double maxZ;
	final double maxX;
	final double maxY;
	final String title;
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
		Optional<CoordinateVector> maZ =
			values.keySet().stream().parallel().max(Comparator.comparingDouble(CoordinateVector::z));
		Optional<CoordinateVector> miZ =
			values.keySet().stream().parallel().min(Comparator.comparingDouble(CoordinateVector::z));
		maxValue = maV.orElse(0.0);
		minValue = miV.orElse(0.0);
		maxX = maX.orElse(new CoordinateVector(2)).x();
		minX = miX.orElse(new CoordinateVector(2)).x();
		maxY = maY.orElse(new CoordinateVector(2)).y();
		minY = miY.orElse(new CoordinateVector(2)).y();
		maxZ = maZ.orElse(new CoordinateVector(2)).z();
		minZ = miZ.orElse(new CoordinateVector(2)).z();
		this.title = title;
	}
	public List<CoordinateVector> getClosestPoints(double slider)
	{
		double currentZ = minZ + slider*(maxZ - minZ);
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
		double maxValueDifference = maxValue - minValue;
		
		for(CoordinateVector co:getClosestPoints(slider))
		{
			int red = (int)((values.get(co))/maxValueDifference*255);
			int green = 0;
			int blue = (int)((maxValue - values.get(co))/maxValueDifference*255);
			Color c = new Color(red,green,blue);
			int x = (int)((co.x() - minX)/(maxX-minX)*(width-150)+75 - pixelWidth/2);
			int y =
				(int)(height - ((co.y() - minY)/(maxY-minY)*(height-150)+75 + pixelHeight/2));
			g.setColor(c);
			g.fillRect(x,y,pixelWidth,pixelHeight);
		}
	}
	
	@Override
	public String title()
	{
		return "3D Plot " + title;
	}
}
