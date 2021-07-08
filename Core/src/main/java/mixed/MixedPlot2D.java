package mixed;

import basic.Plot;
import basic.ScalarFunction;
import basic.ScalarPlot2D;
import linalg.CoordinateVector;

import java.awt.*;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;

public class MixedPlot2D extends ScalarPlot2D
{
	final Map<CoordinateVector, CoordinateVector> velocities;
	double maxV;
	public MixedPlot2D(MixedFunction function, List<CoordinateVector> points, int pointsPerDimension)
	{
		super(function.getPressureFunction(), points, pointsPerDimension);
		velocities = function.getVelocityFunction().valuesInPoints(points);
		OptionalDouble maxVelocity =
			velocities.values().stream().parallel().mapToDouble(CoordinateVector::euclidianNorm).max();
		maxV = maxVelocity.orElse(1);
	}
	
	@Override
	public void drawValues(Graphics g, int width, int height, double slider)
	{
		super.drawValues(g, width, height, slider);
		
		for(Map.Entry<CoordinateVector, CoordinateVector> entry:velocities.entrySet())
		{
			int x = (int)((entry.getKey().x() - minX)/(maxX-minX)*(width-150)+75 );
			int y =
				(int)(height - ((entry.getKey().y() - minY)/(maxY-minY)*(height-150)+75));
			g.setColor(Color.BLACK);
			g.fillOval(x-2,y-2,4,4);
			int vx = (int)(entry.getValue().x()/maxValueDifference*pixelWidth);
			int vy = (int)(entry.getValue().y()/maxValueDifference*pixelHeight);
			g.drawLine(x,y,x+vx,y+vy);
		}
		
	}
	
	@Override
	public String title()
	{
		return "Mixed 2D plot";
	}
	
}
