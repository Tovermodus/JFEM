package mixed;

import basic.ScalarFunction;
import basic.ScalarPlot2D;
import linalg.CoordinateVector;

import java.awt.*;
import java.util.*;
import java.util.List;

public class MixedPlot2D extends ScalarPlot2D
{
	final Map<CoordinateVector, CoordinateVector> velocities;
	double maxV;
	public MixedPlot2D(MixedFunction function, List<CoordinateVector> points, int pointsPerDimension)
	{
		super(function.hasPressureFunction()?function.getPressureFunction(): ScalarFunction.constantFunction(0), points,
			pointsPerDimension);
		if(function.hasVelocityFunction())
			velocities = function.getVelocityFunction().valuesInPoints(points);
		else
			velocities =
				ScalarFunction.constantFunction(0).makeIsotropicVectorFunction().valuesInPoints(points);
		System.out.println(function.hasPressureFunction() + " " + function.hasVelocityFunction());
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
			g.setColor(Color.BLACK);
			g.fillOval(x-2,y-2,4,4);
			int vx = (int)(entry.getValue().x()/(maxV)*pixelWidth);
			int vy = (int)(entry.getValue().y()/(maxV)*pixelHeight);
			g.drawLine(x,y,x+vx,y-vy);
		}
		
	}
	
	@Override
	public String title()
	{
		return "Mixed 2D Plot " + title + " blue is min with " + minValue+" red is max with " + maxValue +
		"max velocity is" + maxV;
	}
	
}
