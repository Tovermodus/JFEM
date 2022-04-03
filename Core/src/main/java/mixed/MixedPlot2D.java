package mixed;

import basic.*;
import linalg.CoordinateVector;

import java.awt.*;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;

public class MixedPlot2D
	extends ScalarPlot2D
{
	Map<CoordinateVector, CoordinateVector> velocities;
	double maxV;
	
	public MixedPlot2D(final MixedFunction function,
	                   final List<CoordinateVector> points,
	                   final int pointsPerDimension)
	{
		super(function.hasPressureFunction() ? function.getPressureFunction() :
		      ScalarFunction.constantFunction(0), points,
		      pointsPerDimension);
		saveVelocities(function, points);
		final OptionalDouble maxVelocity =
			velocities.values()
			          .stream()
			          .parallel()
			          .mapToDouble(CoordinateVector::euclidianNorm)
			          .max();
		maxV = maxVelocity.orElse(1);
	}
	
	public <CT extends Cell<CT, FT>, FT extends Face<CT, FT>> MixedPlot2D(final MixedFunctionOnCells<CT, FT> function,
	                                                                      final FESpace<CT, FT, ?> space,
	                                                                      final int points)
	{
		super(function.hasPressureFunction() ? function.getPressureFunction() :
		      ScalarFunction.constantFunction(0), space.generatePlotPoints(points),
		      points);
		saveVelocities(function, space, points);
		final OptionalDouble maxVelocity =
			velocities.values()
			          .stream()
			          .parallel()
			          .mapToDouble(CoordinateVector::euclidianNorm)
			          .max();
		maxV = maxVelocity.orElse(1);
	}
	
	public MixedPlot2D(final MixedFunction function,
	                   final List<CoordinateVector> points,
	                   final int pointsPerDimension,
	                   final String title)
	{
		super(function.hasPressureFunction() ? function.getPressureFunction() :
		      ScalarFunction.constantFunction(0), points,
		      pointsPerDimension, title);
		saveVelocities(function, points);
		final OptionalDouble maxVelocity =
			velocities.values()
			          .stream()
			          .parallel()
			          .mapToDouble(CoordinateVector::euclidianNorm)
			          .max();
		maxV = maxVelocity.orElse(1);
	}
	
	public void saveVelocities(final MixedFunction function, final List<CoordinateVector> points)
	{
		if (function.hasVelocityFunction())
			velocities = function.getVelocityFunction()
			                     .valuesInPoints(points);
		else
			velocities =
				ScalarFunction.constantFunction(0)
				              .makeIsotropicVectorFunction()
				              .valuesInPoints(points);
	}
	
	public <CT extends Cell<CT, FT>, FT extends Face<CT, FT>> void saveVelocities(final MixedFunctionOnCells<CT, FT> function,
	                                                                              final FESpace<CT, FT, ?> space,
	                                                                              final int points)
	{
		final var poi = space.generatePlotPoints(points);
		if (function.hasVelocityFunction())
			velocities = function.getVelocityFunction()
			                     .valuesInPointsOnCells(poi, space.getCells());
		else
			velocities =
				ScalarFunction.constantFunction(0)
				              .makeIsotropicVectorFunction()
				              .valuesInPoints(poi);
	}
	
	@Override
	public void drawValues(final Graphics g, final int width, final int height, final double slider)
	{
		super.drawValues(g, width, height, slider);
		for (final Map.Entry<CoordinateVector, CoordinateVector> entry : velocities.entrySet())
		{
			final int x = (int) ((entry.getKey()
			                           .x() - min.x()) / (max.x() - min.x()) * (width - 150) + 75);
			final int y = (int) (height - ((entry.getKey()
			                                     .y() - min.y()) / (max.y() - min.y()) * (height - 150) + 75));
			g.setColor(Color.BLACK);
			g.fillOval(x - 2, y - 2, 4, 4);
			final CoordinateVector dir = entry.getValue()
			                                  .mul(1. / entry.getValue()
			                                                 .euclidianNorm());
			final int vx = (int) (dir.x() * pixelWidth);
			final int vy = (int) (dir.y() * pixelHeight);
			g.drawLine(x, y, x + vx, y - vy);
		}
	}
	
	@Override
	public String title()
	{
		return "Mixed 2D Plot " + title + " blue is min with " + minValue + " red is max with " + maxValue +
			"max velocity is" + maxV;
	}
}
