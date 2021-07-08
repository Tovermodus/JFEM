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
	
	@Override
	public void drawValues(Graphics g, int width, int height, double slider)
	{
		super.drawValues(g, width, height, slider);
		for(Map.Entry<CoordinateVector, CoordinateVector> entry:velocities.entrySet())
		{
			int x = (int)((entry.getKey().x() - min.x())/(max.x()-min.x())*(width-150)+75 );
			int y =	(int)(height - ((entry.getKey().y() - min.y())/(max.y()-min.y())*(height-150)+75));
			int green = (int)(entry.getValue().z()/(maxValue - minValue)*255);
			if(green>255)
				green = 255;
			g.setColor(new Color(0,green,0));
			g.fillOval(x-2,y-2,4,4);
			int vx = (int)(entry.getValue().x()/(maxValue - minValue)*pixelWidth);
			int vy = (int)(entry.getValue().y()/(maxValue - minValue)*pixelHeight);
			g.drawLine(x,y,x+vx,y+vy);
		}
		
	}
	
	@Override
	public String title()
	{
		return "Mixed 2D plot";
	}
	
}
