package mixed;

import linalg.CoordinateVector;

import java.util.Map;
import java.util.OptionalDouble;

public class MixedPlot2DTime extends MixedPlot3D
{
	public MixedPlot2DTime(final Map<CoordinateVector, Double> pressures, final Map<CoordinateVector, CoordinateVector> velocities, final int pointsPerDimension)
	{
		this(pressures, velocities, pointsPerDimension, "unnamed");
	}
	
	public MixedPlot2DTime(final Map<CoordinateVector, Double> pressures, final Map<CoordinateVector, CoordinateVector> velocities, final int pointsPerDimension, final String title)
	{
		super(pressures, velocities, pointsPerDimension, title);
		final OptionalDouble maxVelocity = velocities
			.values()
			.stream()
			.parallel()
			.mapToDouble(c -> Math.sqrt(c.x() * c.x() + c.y() * c.y()))
			.max();
		maxV = maxVelocity.orElse(1);
	}
	
	@Override
	public String title()
	{
		return "Mixed 2D Plot in Time " + title + " blue is min with " + minValue + " red is max with " + origMaxVal + "\nmax velocity is" + maxV + "time is " + z;
	}
}
