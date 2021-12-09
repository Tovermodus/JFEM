package distorted.geometry;

import basic.DoubleCompare;
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

public class RingGrid
	implements DistortedGrid
{
	public final int dimension;
	private final CoordinateVector centerPoint;
	private final double innerRadius;
	private final double outerRadius;
	public ImmutableSet<DistortedCell> cells;
	public ImmutableSet<DistortedFace> faces;
	
	public RingGrid(final CoordinateVector centerPoint, final double innerRadius, final double outerRadius,
	                final int refinements)
	{
		this.dimension = centerPoint.getLength();
		this.centerPoint = centerPoint;
		this.innerRadius = innerRadius;
		this.outerRadius = outerRadius;
		
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
			final double innerCircumference = 2 * Math.PI * innerRadius;
			final int nCells = (int) (0.8 * innerCircumference / (outerRadius - innerRadius));
			for (int i = 0; i < nCells; i++)
			{
				final double firstDegree = Math.PI * 2 / nCells * i;
				final double secondDegree = Math.PI * 2 / nCells * (i + 1);
				final CoordinateVector inner_first
					= CoordinateVector.fromValues(Math.cos(firstDegree),
					                              Math.sin(firstDegree))
					                  .mul(innerRadius)
					                  .add(centerPoint);
				final CoordinateVector inner_second
					= CoordinateVector.fromValues(Math.cos(secondDegree),
					                              Math.sin(secondDegree))
					                  .mul(innerRadius)
					                  .add(centerPoint);
				final CoordinateVector outer_first
					= CoordinateVector.fromValues(Math.cos(firstDegree),
					                              Math.sin(firstDegree))
					                  .mul(outerRadius)
					                  .add(centerPoint);
				final CoordinateVector outer_second
					= CoordinateVector.fromValues(Math.cos(secondDegree),
					                              Math.sin(secondDegree))
					                  .mul(outerRadius)
					                  .add(centerPoint);
				final List<CoordinateVector> vertices = List.of(inner_first, outer_first, outer_second,
				                                                inner_second);
				addCell(vertices, ret);
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
			final boolean isOnInnerBoundary = isOnInnerBoundary(centerPoint, vertices);
			final boolean isOnOuterBoundary = isOnOuterBoundary(centerPoint, vertices);
			final DistortedFace face = new DistortedFace(vertices, isOnInnerBoundary || isOnOuterBoundary);
			for (final DistortedCell c : DistortedGrid.getCellsBelongingToFace(face, neighbouringCells))
			{
				c.faces.add(face);
				face.addCell(c);
			}
			faces.add(face);
		}
	}
	
	private boolean isOnInnerBoundary(final CoordinateVector centerPoint, final CoordinateVector[] vertices)
	{
		return Arrays.stream(vertices)
		             .map(centerPoint::sub)
		             .mapToDouble(CoordinateVector::euclidianNorm)
		             .allMatch(norm -> DoubleCompare.almostEqual(norm,
		                                                         innerRadius));
	}
	
	private boolean isOnOuterBoundary(final CoordinateVector centerPoint, final CoordinateVector[] vertices)
	{
		return Arrays.stream(vertices)
		             .map(centerPoint::sub)
		             .mapToDouble(CoordinateVector::euclidianNorm)
		             .allMatch(norm -> DoubleCompare.almostEqual(norm,
		                                                         outerRadius));
	}
	
	@Override
	public List<DistortedCell> partitionCell(final DistortedCell cell)
	{
		final List<DistortedCell> ret = new ArrayList<>();
		final CoordinateVector[] vertices = cell.vertices;
		if (dimension == 2)
		{
			final DistortedFace bottom = cell.getFaceFromVertexNumbers(0, 1);
			CoordinateVector bottomCenter = bottom.center();
			final DistortedFace right = cell.getFaceFromVertexNumbers(1, 2);
			CoordinateVector rightCenter = right.center();
			final DistortedFace top = cell.getFaceFromVertexNumbers(2, 3);
			CoordinateVector topCenter = top.center();
			final DistortedFace left = cell.getFaceFromVertexNumbers(3, 0);
			CoordinateVector leftCenter = left.center();
			if (bottom.isBoundaryFace())
				bottomCenter = centerPoint.add(bottomCenter.sub(centerPoint)
				                                           .normalize()
				                                           .mul(bottom.getVertices()
				                                                      .get(0)
				                                                      .euclidianNorm()));
			if (right.isBoundaryFace())
				rightCenter = centerPoint.add(rightCenter.sub(centerPoint)
				                                         .normalize()
				                                         .mul(right.getVertices()
				                                                   .get(0)
				                                                   .euclidianNorm()));
			if (top.isBoundaryFace())
				topCenter = centerPoint.add(topCenter.sub(centerPoint)
				                                     .normalize()
				                                     .mul(top.getVertices()
				                                             .get(0)
				                                             .euclidianNorm()));
			if (left.isBoundaryFace())
				leftCenter = centerPoint.add(leftCenter.sub(centerPoint)
				                                       .normalize()
				                                       .mul(left.getVertices()
				                                                .get(0)
				                                                .euclidianNorm()));
			final CoordinateVector center = mean(bottomCenter, rightCenter, topCenter, leftCenter);
			System.out.println();
			System.out.println("ll " + vertices[0]);
			System.out.println("lr " + vertices[1]);
			System.out.println("ur " + vertices[2]);
			System.out.println("ul " + vertices[3]);
			System.out.println("bc " + bottomCenter);
			System.out.println("rc " + rightCenter);
			System.out.println("tc " + topCenter);
			System.out.println("lc " + leftCenter);
			System.out.println("cc " + center);
			
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
					                for (int i = 0; i < dimension; i++)
					                {
						                coords.set(1. / (pointsPerCell - 1) *
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
		if (getDimension() != 2)
			throw new UnsupportedOperationException("Not yetimplemented");
		final List<CoordinateVector> ret = new ArrayList<>(resolution * resolution);
		for (int i = 0; i < resolution; i++)
		{
			for (int j = 0; j < resolution; j++)
			{
				final double ang = 1.0 * i / resolution * Math.PI * 2;
				final double rad = innerRadius + 1.0 * (outerRadius - innerRadius) * j / resolution;
				final CoordinateVector out =
					CoordinateVector.fromValues(Math.cos(ang), Math.sin(ang))
					                .mul(rad)
					                .add(centerPoint);
				if (cells.stream()
				         .anyMatch(cell -> cell.isInCell(out)))
					ret.add(new CoordinateVector(centerPoint.add(out)));
			}
		}
		return ret;
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
					                for (int i = 0; i < dimension; i++)
					                {
						                coords.set(1. / (pointsPerCell - 1) *
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
