package mixed;

import basic.ScalarPlot2D;
import basic.ScalarPlot3D;
import linalg.CoordinateVector;

import java.awt.*;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;

public class MixedPlot3D extends ScalarPlot3D
{
	final Map<CoordinateVector, CoordinateVector> velocities;
	double maxV;
	public MixedPlot3D(MixedFunction function, List<CoordinateVector> points, int pointsPerDimension)
	{
		super(function.getPressureFunction(), points, pointsPerDimension);
		velocities = function.getVelocityFunction().valuesInPoints(points);
		OptionalDouble maxVelocity =
			velocities.values().stream().parallel().mapToDouble(CoordinateVector::euclidianNorm).max();
		maxV = maxVelocity.orElse(1);
	}
	public MixedPlot3D(Map<CoordinateVector, Double> pressures,
	                   Map<CoordinateVector, CoordinateVector> velocities, int pointsPerDimension)
	{
		super(pressures, pointsPerDimension, "");
		this.velocities = velocities;
		OptionalDouble maxVelocity =
			velocities.values().stream().parallel().mapToDouble(CoordinateVector::euclidianNorm).max();
		maxV = maxVelocity.orElse(1);
	}
	@Override
	public void drawValues(Graphics g, int width, int height, double slider)
	{
		super.drawValues(g, width, height, slider);
		double velocityScaling = Math.max(Math.sqrt((maxValue - minValue)*maxV), maxV);
		for(CoordinateVector c:getClosestPoints(slider))
		{
			int x = (int)((c.x() - min.x())/(max.x()-min.x())*(width-150)+75 );
			int y =	(int)(height - ((c.y() - min.y())/(max.y()-min.y())*(height-150)+75));
			int green = (int)((velocities.get(c).z() - min.z())/velocityScaling*255);
			if(green>255)
				green = 255;
			if(green<0)
				green = 0;
			g.setColor(new Color(0,green,0));
			g.fillOval(x-2,y-2,4,4);
			int vx = (int)(velocities.get(c).x() /velocityScaling*pixelWidth);
			int vy = (int)(velocities.get(c).y() /velocityScaling*pixelHeight);
			g.drawLine(x,y,x+vx,y-vy);
		}
		
	}
	
	@Override
	public String title()
	{
		return "Mixed 2D plot";
	}
	
}
