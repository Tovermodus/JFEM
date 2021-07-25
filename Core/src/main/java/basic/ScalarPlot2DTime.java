package basic;

import com.google.common.collect.Iterables;
import linalg.CoordinateVector;

import java.awt.*;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Optional;

public class ScalarPlot2DTime extends ScalarPlot3D
{
	public ScalarPlot2DTime(Map<CoordinateVector, Double> values, int pointsPerDimension,
	                        String title)
	{
		super(values, pointsPerDimension, title);
	}
	@Override
	public String title()
	{
		return "2D Plot with time " + title + " blue is min with " + minValue+" red is max with " + maxValue;
	}
}
