package mixed;

import basic.ScalarPlot2D;
import linalg.CoordinateVector;

import java.awt.*;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;

public class MixedPlot2DTime extends MixedPlot3D
{
	public MixedPlot2DTime(Map<CoordinateVector, Double> pressures,
	                       Map<CoordinateVector, CoordinateVector> velocities,
	                       int pointsPerDimension)
	{
		super(pressures, velocities, pointsPerDimension);
	}
	
}
