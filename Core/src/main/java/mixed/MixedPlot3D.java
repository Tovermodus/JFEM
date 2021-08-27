package mixed;

import basic.ScalarFunction;
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
	
	public MixedPlot3D(final MixedFunction function, final List<CoordinateVector> points, final int pointsPerDimension)
	{
		super(function.hasPressureFunction() ? function.getPressureFunction() :
		      ScalarFunction.constantFunction(0), points, pointsPerDimension);
		if (function.hasVelocityFunction()) velocities = function.getVelocityFunction().valuesInPoints(points);
		else velocities = ScalarFunction
			.constantFunction(0)
			.makeIsotropicVectorFunction()
			.valuesInPoints(points);
		final OptionalDouble maxVelocity = velocities
			.values()
			.stream()
			.parallel()
			.mapToDouble(CoordinateVector::euclidianNorm)
			.max();
		maxV = maxVelocity.orElse(1);
	}
	
	public MixedPlot3D(final Map<CoordinateVector, Double> pressures, final Map<CoordinateVector, CoordinateVector> velocities, final int pointsPerDimension, final String title)
	{
		super(pressures, pointsPerDimension, title);
		this.velocities = velocities;
		final OptionalDouble maxVelocity = velocities
			.values()
			.stream()
			.parallel()
			.mapToDouble(CoordinateVector::euclidianNorm)
			.max();
		maxV = maxVelocity.orElse(1);
	}
	
	@Override
	public void drawValues(final Graphics g, final int width, final int height, final double slider)
	{
		super.drawValues(g, width, height, slider);
		for (final CoordinateVector c : getClosestPoints(slider))
		{
			final int x = (int) ((c.x() - min.x()) / (max.x() - min.x()) * (width - 150) + 75);
			final int y = (int) (height - ((c.y() - min.y()) / (max.y() - min.y()) * (height - 150) + 75));
			g.setColor(Color.BLACK);
			g.fillOval(x - 2, y - 2, 4, 4);
			final int vx = (int) (velocities.get(c).x() / (maxV) * pixelWidth);
			final int vy = (int) (velocities.get(c).y() / (maxV) * pixelHeight);
			g.drawLine(x, y, x + vx, y - vy);
		}
	}
	
	@Override
	public String title()
	{
		
		return "Mixed 3D Plot " + title + " blue is min with " + minValue + " red is max with " + maxValue + "\nmax velocity is" + maxV + " z is " + z;
	}
}
