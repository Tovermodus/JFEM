package tensorproduct.geometry;

import basic.PerformanceArguments;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

public class CartesianGrid
{
	public final CoordinateVector startCoordinates;
	public final CoordinateVector endCoordinates;
	public final IntCoordinates cellsPerDimension;
	public ImmutableList<ImmutableList<Cell1D>> cells1D;
	public final HashMap<IntCoordinates, TPCell> cellsByCoordinates;
	public ImmutableList<TPCell> cells;
	public ImmutableList<TPFace> faces;
	public final int dimension;
	
	private void createCells1D()
	{
		ArrayList<ArrayList<Cell1D>> mutableCells1D = new ArrayList<>(3);
		for (int i = 0; i < dimension; i++)
		{
			double end = endCoordinates.at(i);
			double start = startCoordinates.at(i);
			int count = cellsPerDimension.get(i);
			double width = (end - start) / count;
			ArrayList<Cell1D> cellsThisDimension = new ArrayList<>(count);
			for (int j = 0; j < count; j++)
			{
				cellsThisDimension.add(new Cell1D(j * width + start, (j + 1) * width + start));
			}
			mutableCells1D.add(cellsThisDimension);
		}
		ArrayList<ImmutableList<Cell1D>> intermediateCopy = new ArrayList<>();
		for (ArrayList<Cell1D> l : mutableCells1D)
			intermediateCopy.add(ImmutableList.copyOf(l));
		cells1D = ImmutableList.copyOf(intermediateCopy);
		
	}
	
	public CartesianGrid(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                     IntCoordinates cellsPerDimension)
	{
		this.startCoordinates = startCoordinates;
		this.endCoordinates = endCoordinates;
		this.cellsPerDimension = cellsPerDimension;
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (startCoordinates.getLength() != endCoordinates.getLength() || endCoordinates.getLength() != cellsPerDimension.getDimension())
				throw new IllegalArgumentException("Dimensions of coordinates or cell counts do not match");
		}
		dimension = startCoordinates.getLength();
		createCells1D();
		cellsByCoordinates = new HashMap<>(cellsPerDimension.size());
		createCells();
		createFaces();
	}
	
	private void createFaces()
	{
		ArrayList<TPFace> mutableFaces = new ArrayList<>(dimension * cellsPerDimension.size());
		for (int normalDirection = 0; normalDirection < dimension; normalDirection++)
		{
			IntCoordinates numberOfFaces = new IntCoordinates(cellsPerDimension);
			for (int i = 0; i < dimension; i++)
				if (i == normalDirection)
					numberOfFaces = numberOfFaces.get_modified(i, cellsPerDimension.get(i) + 1);
			for (IntCoordinates c : numberOfFaces.range())
			{
				TPFace newFace = getNewFace(normalDirection, c);
				combineFaceWithNeighbouringCells(normalDirection, c, newFace);
				mutableFaces.add(newFace);
			}
		}
		faces = ImmutableList.copyOf(mutableFaces);
	}
	
	@NotNull
	private TPFace getNewFace(int normalDirection, IntCoordinates c)
	{
		List<Cell1D> componentCells = new ArrayList<>(dimension - 1);
		for (int i = 0; i < dimension; i++)
		{
			if (i != normalDirection)
				componentCells.add(cells1D.get(i).get(c.get(i)));
		}
		boolean boundaryFace = false;
		double otherCoordinate = 0;
		int layer = c.get(normalDirection);
		int topLayer = cellsPerDimension.get(normalDirection);
		if (layer == 0)
			boundaryFace = true;
		if (layer < topLayer)
			otherCoordinate =
				cells1D.get(normalDirection).get(layer).getStart();
		if (layer == topLayer)
		{
			boundaryFace = true;
			otherCoordinate = cells1D.get(normalDirection).get(layer - 1).getEnd();
		}
		return new TPFace(componentCells, normalDirection, otherCoordinate,
			boundaryFace);
	}
	
	
	private void combineFaceWithNeighbouringCells(int normalDirection, IntCoordinates c, TPFace newFace)
	{
		int topLayer = cellsPerDimension.get(normalDirection);
		int layer = c.get(normalDirection);
		if (layer > 0)
		{
			TPCell lowerCell = cellsByCoordinates.get(c.get_modified(normalDirection, layer - 1));
			newFace.cells.add(lowerCell);
			lowerCell.faces.add(newFace);
		}
		if (layer < topLayer)
		{
			TPCell upperCell = cellsByCoordinates.get(c.get_modified(normalDirection, layer));
			newFace.cells.add(upperCell);
			upperCell.faces.add(newFace);
		}
	}
	
	private void createCells()
	{
		ArrayList<TPCell> mutableCells = new ArrayList<>(cellsPerDimension.size());
		for (IntCoordinates c : cellsPerDimension.range())
		{
			List<Cell1D> componentCells = new ArrayList<>(dimension);
			for (int i = 0; i < dimension; i++)
				componentCells.add(cells1D.get(i).get(c.get(i)));
			TPCell newCell = new TPCell(componentCells);
			cellsByCoordinates.put(c, newCell);
			mutableCells.add(newCell);
		}
		cells = ImmutableList.copyOf(mutableCells);
	}
	
	public List<CoordinateVector> generatePlotPoints(int resolution)
	{
		IntCoordinates pointsPerDimension = IntCoordinates.repeat(resolution, dimension);
		return pointsPerDimension.range().stream().map(c ->
		{
			CoordinateVector coords = new CoordinateVector(dimension);
			for (int i = 0; i < dimension; i++)
			{
				coords.set(startCoordinates.at(i) + (endCoordinates.at(i) - startCoordinates.at(i)) / (resolution-1) * c.get(i), i);
			}
			return coords;
		}).collect(Collectors.toList());
	}
}
