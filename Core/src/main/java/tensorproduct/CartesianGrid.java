package tensorproduct;

import com.google.common.collect.ImmutableList;
import linalg.CoordinateVector;

import java.util.List;

public class CartesianGrid
{
	ImmutableList<ImmutableList<Cell1D>> cells1D;
	ImmutableList<TPCell> cells;
	ImmutableList<TPFace> faces;
	
	public CartesianGrid(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                         List<Integer> cellsPerDimension)
	{
	
	}
}
