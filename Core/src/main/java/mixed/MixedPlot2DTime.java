package mixed;

import linalg.CoordinateVector;

import java.util.Map;
import java.util.OptionalDouble;

public class MixedPlot2DTime extends MixedPlot3D
{
	public MixedPlot2DTime(Map<CoordinateVector, Double> pressures,
	                       Map<CoordinateVector, CoordinateVector> velocities,
	                       int pointsPerDimension)
	{
		super(pressures, velocities, pointsPerDimension);
		OptionalDouble maxVelocity =
			velocities.values().stream().parallel().mapToDouble(c->Math.sqrt(c.x()*c.x() + c.y()*c.y())).max();
		maxV = maxVelocity.orElse(1);
	}
	
	@Override
	public String title()
	{
		return "Mixed 2D Plot in Time " + title + " blue is min with " + minValue+" red is max with " + maxValue +
			"max velocity is" + maxV;
	}
}
