package tensorproduct;

import basic.PerformanceArguments;
import com.google.common.collect.ImmutableList;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.checkerframework.checker.units.qual.A;

import java.util.ArrayList;
import java.util.List;

public class CartesianGrid
{
	private final CoordinateVector startCoordinates;
	private final CoordinateVector endCoordinates;
	private final IntCoordinates cellsPerDimension;
	ImmutableList<ImmutableList<Cell1D>> cells1D;
	ImmutableList<TPCell> cells;
	ImmutableList<TPFace> faces;
	int dimension;
	private void createCells1D()
	{
		ArrayList<ArrayList<Cell1D>> mutableCells1D = new ArrayList<>(3);
		for(int i = 0; i < dimension; i++)
		{
			double end = endCoordinates.at(i);
			double start = startCoordinates.at(i);
			int count = cellsPerDimension.get(i);
			double width = (end - start)/count;
			ArrayList<Cell1D> cellsThisDimension = new ArrayList<>(count);
			for (int j = 0; j < count; j++)
			{
				cellsThisDimension.add(new Cell1D(j*width, (j+1)*width));
			}
			mutableCells1D.add(cellsThisDimension);
		}
		ArrayList<ImmutableList<Cell1D>> intermediateCopy = new ArrayList<>();
		for(ArrayList<Cell1D> l: mutableCells1D)
			intermediateCopy.add(ImmutableList.copyOf(l));
		cells1D = ImmutableList.copyOf(intermediateCopy);
		
	}
	public CartesianGrid(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                         IntCoordinates cellsPerDimension)
	{
		this.startCoordinates = startCoordinates;
		this.endCoordinates = endCoordinates;
		this.cellsPerDimension = cellsPerDimension;
		if(PerformanceArguments.getInstance().executeChecks)
		{
			if(startCoordinates.getLength() != endCoordinates.getLength() || endCoordinates.getLength() != cellsPerDimension.size())
				throw new IllegalArgumentException("Dimensions of coordinates or cell counts do not " +
					"match");
		}
		dimension = startCoordinates.getLength();
		createCells1D();
	}
}
