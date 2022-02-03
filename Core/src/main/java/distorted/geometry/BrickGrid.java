package distorted.geometry;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableCollection;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import io.vavr.Tuple2;
import kotlin.Pair;
import linalg.CoordinateVector;
import linalg.IntCoordinates;

import java.util.*;
import java.util.stream.Collectors;

public class BrickGrid
	implements DistortedGrid
{
	public final int dimension;
	private final CoordinateVector centerPoint;
	private final double width;
	private final double height;
	public ImmutableSet<DistortedCell> cells;
	public ImmutableSet<DistortedFace> faces;
	
	public BrickGrid(final CoordinateVector centerPoint, final double width, final double height,
	                 final int refinements)
	{
		this.dimension = centerPoint.getLength();
		this.centerPoint = centerPoint;
		this.width = width;
		this.height = height;
		
		final Pair<ImmutableSet<DistortedCell>, ImmutableSet<DistortedFace>> grid = createGrid(refinements);
		cells = grid.component1();
		faces = grid.component2();
	}
	
	private static void addCell(final List<DistortedCell> list, final CoordinateVector... verticesArray)
	{
		final List<CoordinateVector> vertices = List.of(verticesArray);
		final CoordinateVector rating;
		rating = CoordinateVector.fromValues(1, 1);
		final Optional<CoordinateVector> minRating = vertices.stream()
		                                                     .min(Comparator.comparingDouble(rating::inner));
		if (minRating.isEmpty())
			throw new IllegalStateException("called with no vertices");
		final int minIndex = vertices.indexOf(minRating.get());
		final List<CoordinateVector> orderedVertices = new ArrayList<>(4);
		for (int j = 0; j < 4; j++)
		{
			orderedVertices.add(vertices.get((minIndex + j) % 4));
		}
		final DistortedCell cell =
			new DistortedCell(orderedVertices.toArray(CoordinateVector[]::new));
		list.add(cell);
	}
	
	private static void addCell(final List<CoordinateVector> vertices,
	                            final Multimap<DistortedCell, DistortedCell> map)
	{
		final CoordinateVector rating;
		rating = CoordinateVector.fromValues(1, 1);
		final Optional<CoordinateVector> minRating = vertices.stream()
		                                                     .min(Comparator.comparingDouble(rating::inner));
		if (minRating.isEmpty())
			throw new IllegalStateException("called with no vertices");
		final int minIndex = vertices.indexOf(minRating.get());
		final List<CoordinateVector> orderedVertices = new ArrayList<>(4);
		for (int j = 0; j < 4; j++)
		{
			orderedVertices.add(vertices.get((minIndex + j) % 4));
		}
		final DistortedCell cell =
			new DistortedCell(orderedVertices.toArray(CoordinateVector[]::new));
		map.put(cell, cell);
	}
	
	@Override
	public Multimap<DistortedCell, DistortedCell> generateCells()
	{
		final Multimap<DistortedCell, DistortedCell> ret = HashMultimap.create();
		if (dimension == 2)
		{
			final double minSide = Math.min(width, height);
			final double cellWidthEst = minSide / 1.9999;
			final int xCells = (int) (width / cellWidthEst) + 1;
			final int yCells = (int) (height / cellWidthEst) + 1;
			final double cellWidth = width / xCells;
			final double cellHeight = height / yCells;
			final double xStart = centerPoint.x() - width / 2;
			final double yStart = centerPoint.y() - height / 2;
			for (int i = 0; i < xCells; i++)
			{
				for (int j = 0; j < yCells; j++)
				{
					final CoordinateVector ll = CoordinateVector.fromValues(xStart + i * cellWidth,
					                                                        yStart + j * cellHeight);
					final CoordinateVector lr
						= CoordinateVector.fromValues(xStart + (i + 1) * cellWidth,
						                              yStart + j * cellHeight);
					final CoordinateVector ul = CoordinateVector.fromValues(xStart + i * cellWidth,
					                                                        yStart + (j + 1) * cellHeight);
					final CoordinateVector ur
						= CoordinateVector.fromValues(xStart + (i + 1) * cellWidth,
						                              yStart + (j + 1) * cellHeight);
					final List<CoordinateVector> vertices = List.of(ll,
					                                                lr,
					                                                ur,
					                                                ul);
					System.out.println(i + " " + j + " " + vertices);
					addCell(vertices, ret);
				}
			}
			return ret;
		} else throw new IllegalArgumentException("Only supports dimension 2");
	}
	
	@Override
	public void createFaces(final HashSet<DistortedFace> faces, final DistortedCell cell,
	                        final Collection<DistortedCell> neighbouringCells)
	{
		for (final CoordinateVector[] vertices : cell.getSides())
		{
			final boolean isOnBound = isOnBoundary(vertices);
			final DistortedFace face = new DistortedFace(vertices, isOnBound);
			for (final DistortedCell c : DistortedGrid.getCellsBelongingToFace(face, neighbouringCells))
			{
				c.faces.add(face);
				face.addCell(c);
			}
			faces.add(face);
		}
	}
	
	private boolean isOnBoundary(final CoordinateVector[] vertices)
	{
		for (final CoordinateVector vertex : vertices)
		{
			if (vertex.x() == centerPoint.x() - width / 2)
				return true;
			if (vertex.x() == centerPoint.x() + width / 2)
				return true;
			if (vertex.y() == centerPoint.y() - height / 2)
				return true;
			if (vertex.y() == centerPoint.y() + height / 2)
				return true;
		}
		return false;
	}
	
	@Override
	public List<DistortedCell> partitionCell(final DistortedCell cell)
	{
		final List<DistortedCell> ret = new ArrayList<>();
		final CoordinateVector[] vertices = cell.vertices;
		if (dimension == 2)
		{
			final DistortedFace bottom = cell.getFaceFromVertexNumbers(0, 1);
			final CoordinateVector bottomCenter = bottom.center();
			final DistortedFace right = cell.getFaceFromVertexNumbers(1, 2);
			final CoordinateVector rightCenter = right.center();
			final DistortedFace top = cell.getFaceFromVertexNumbers(2, 3);
			final CoordinateVector topCenter = top.center();
			final DistortedFace left = cell.getFaceFromVertexNumbers(3, 0);
			final CoordinateVector leftCenter = left.center();
			final CoordinateVector center = mean(bottomCenter, rightCenter, topCenter, leftCenter);
			addCell(ret, vertices[0], bottomCenter, center, leftCenter);
			addCell(ret, bottomCenter, vertices[1], rightCenter, center);
			addCell(ret, leftCenter, center, topCenter, vertices[3]);
			addCell(ret, center, rightCenter, vertices[2], topCenter);
			return ret;
		}
		if (dimension == 3)
		{
			throw new UnsupportedOperationException("Not Yet implemented");
		}
		throw new IllegalStateException("Dimension Wrong");
	}
	
	@Override
	public List<CoordinateVector> generatePlotPoints(final int resolution)
	{
		final int pointsPerCell = resolution / cells.size();
		final List<CoordinateVector> ret = new ArrayList<>();
		final IntCoordinates pointsPerDimension = IntCoordinates.repeat(pointsPerCell, dimension);
		for (final DistortedCell cell : cells)
		{
			ret.addAll(pointsPerDimension
				           .range()
				           .stream()
				           .map(c ->
				                {
					                final CoordinateVector coords = new CoordinateVector(dimension);
					                final double offset = 1. / (3 * pointsPerCell);
					                final double space = (1. - 2 * offset) / (pointsPerCell - 1);
					                for (int i = 0; i < dimension; i++)
					                {
						                coords.set(offset + space *
							                (c.get(i)), i);
					                }
					                return coords;
				                })
				           .map(cell::transformFromReferenceCell)
				           .collect(Collectors.toList()));
		}
		return ret;
	}
	
	@Override
	public List<CoordinateVector> generateIsotropicPlotPoints(final int resolution)
	{
		return generatePlotPoints(resolution);
	}
	
	@Override
	public List<Tuple2<DistortedCell, CoordinateVector>> generateReferencePlotPoints(final int resolution)
	{
		final int pointsPerCell = resolution / cells.size();
		final List<Tuple2<DistortedCell, CoordinateVector>> ret = new ArrayList<>();
		final IntCoordinates pointsPerDimension = IntCoordinates.repeat(pointsPerCell, dimension);
		for (final DistortedCell cell : cells)
		{
			ret.addAll(pointsPerDimension
				           .range()
				           .stream()
				           .map(c ->
				                {
					                final CoordinateVector coords = new CoordinateVector(dimension);
					                final double offset = 1. / (3 * pointsPerCell);
					                final double space = (1. - 2 * offset) / (pointsPerCell - 1);
					                for (int i = 0; i < dimension; i++)
					                {
						                coords.set(offset + space *
							                (c.get(i)), i);
					                }
					                return coords;
				                })
				           .map(c -> new Tuple2<>(cell, c))
				           .collect(Collectors.toList()));
		}
		return ret;
	}
	
	@Override
	public ImmutableCollection<DistortedFace> getFaces()
	{
		return faces;
	}
	
	@Override
	public ImmutableCollection<DistortedCell> getCells()
	{
		return cells;
	}
	
	@Override
	public int getDimension()
	{
		return dimension;
	}
}
