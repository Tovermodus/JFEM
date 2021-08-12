package distorted.geometry;

import basic.FaceWithReferenceFace;
import basic.PerformanceArguments;
import basic.VectorFunction;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import linalg.AffineTransformation;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPFace;

import java.util.*;

public class DistortedFace implements FaceWithReferenceFace<DistortedCell, DistortedFace>
{
	final Map<DistortedCell, AffineTransformation> cellsWithTransformation;
	final TPFace referenceFace;
	final int dimension;
	private final CoordinateVector[] vertices;
	private final boolean isBoundaryFace;
	private final VectorFunction normal;
	private final CoordinateVector normalVector;
	
	public DistortedFace(final CoordinateVector[] vertices, final boolean isBoundaryFace)
	{
		this.vertices = vertices;
		dimension = vertices[0].getLength();
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (vertices.length != (int) Math.pow(2, dimension - 1))
				throw new IllegalArgumentException("Wrong number of vertices");
		}
		this.isBoundaryFace = isBoundaryFace;
		cellsWithTransformation = new HashMap<>();
		if (dimension == 2)
		{
			final CoordinateVector tangent = vertices[1].sub(vertices[0]);
			normalVector = CoordinateVector.fromValues(tangent.at(1), -tangent.at(0));
		} else
		{
			normalVector = vertices[1]
				.sub(vertices[0])
				.cross(vertices[2].sub(vertices[1]));
		}
		normalVector.mulInPlace(1. / normalVector.euclidianNorm());
		normal = new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return dimension;
			}
			
			@Override
			public int getDomainDimension()
			{
				return dimension;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				
				return normalVector;
			}
		};
		referenceFace = TPFace.unitHyperCubeFace(dimension, isBoundaryFace);
	}
	
	@Override
	public int getDimension()
	{
		return dimension;
	}
	
	@Override
	public String toString()
	{
		return "DistortedFace{" +
			"vertices=" + Arrays.toString(vertices) +
			", normalvector=" + normalVector +
			", isBoundaryFace=" + isBoundaryFace +
			'}';
	}
	
	@Override
	public ImmutableSet<DistortedCell> getCells()
	{
		return ImmutableSet.copyOf(cellsWithTransformation.keySet());
	}
	
	public boolean isOnCell(final DistortedCell cell)
	{
		//System.out.println(Arrays.toString(Arrays.stream(vertices).mapToInt(cell::getPositionOfVertex)
		// .toArray()));
		return Arrays
			.stream(vertices)
			.mapToInt(cell::getPositionOfVertex)
			.noneMatch(i -> i == -1);
	}
	
	public AffineTransformation getTransformationFromReferenceFaceToFaceOfReferenceCell(final DistortedCell cell)
	{
		if (!isOnCell(cell))
			throw new IllegalArgumentException("Face does not belong to cell");
		final IntCoordinates vertexNumbers =
			new IntCoordinates(Arrays
				                   .stream(vertices)
				                   .mapToInt(cell::getPositionOfVertex)
				                   .toArray());
		final DistortedCell referenceCell = cell.getReferenceCell();
		final CoordinateVector vector = referenceCell.vertices[vertexNumbers.get(0)];
		final CoordinateMatrix matrix = new CoordinateMatrix(dimension, dimension);
		if (dimension == 2)
		{
			final CoordinateVector mapsFromX = referenceCell.vertices[vertexNumbers.get(1)].sub(vector);
			matrix.addColumn(mapsFromX, 0);
			matrix.addColumn(CoordinateVector
				                 .repeat(1, 2)
				                 .sub(mapsFromX.componentWise(Math::abs)), 1);
		} else if (dimension == 3)
		{
			final CoordinateVector mapsFromY = referenceCell.vertices[vertexNumbers.get(1)].sub(vector);
			final CoordinateVector mapsFromX = referenceCell.vertices[vertexNumbers.get(3)].sub(vector);
			matrix.addColumn(mapsFromX, 0);
			matrix.addColumn(mapsFromY, 1);
			matrix.addColumn(CoordinateVector
				                 .repeat(1, 3)
				                 .sub(mapsFromX.componentWise(Math::abs))
				                 .sub(mapsFromY.componentWise(Math::abs)), 2);
		} else
			throw new IllegalStateException("Dimension must be 2 or 3");
		return new AffineTransformation(matrix, vector);
	}
	
	void addCell(final DistortedCell c)
	{
		if (isOnCell(c))
			cellsWithTransformation.put(c, getTransformationFromReferenceFaceToFaceOfReferenceCell(c));
		else throw new IllegalArgumentException("Face is not on cell");
	}
	
	@Override
	public boolean isBoundaryFace()
	{
		return isBoundaryFace;
	}
	
	@Override
	public VectorFunction getNormal()
	{
		return normal;
	}
	
	@Override
	public DistortedCell getNormalDownstreamCell()
	{
		for (final DistortedCell c : getCells())
		{
			if (center()
				.sub(c.center())
				.inner(normalVector) < 0)
				return c;
		}
		return null;
	}
	
	@Override
	public DistortedCell getNormalUpstreamCell()
	{
		for (final DistortedCell c : getCells())
		{
			if (center()
				.sub(c.center())
				.inner(normalVector) > 0)
				return c;
		}
		return null;
	}
	
	@Override
	public CoordinateVector center()
	{
		return Arrays
			.stream(vertices)
			.reduce(new CoordinateVector(dimension), CoordinateVector::add)
			.mul(
				1. / vertices.length);
	}
	
	@Override
	public boolean equals(final Object o)
	{
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		final DistortedFace that = (DistortedFace) o;
		if (isBoundaryFace != that.isBoundaryFace)
			return false;
		for (final CoordinateVector vertex : vertices)
		{
			boolean otherHasVertex = false;
			for (int j = 0; j < vertices.length; j++)
				if (vertex.almostEqualMute(that.vertices[j]))
				{
					otherHasVertex = true;
					break;
				}
			if (!otherHasVertex)
			{
				return false;
			}
		}
		return true;
	}
	
	@Override
	public int hashCode()
	{
		int result = Objects.hash(isBoundaryFace, dimension);
		result = 31 * result + Arrays
			.stream(vertices)
			.mapToInt(CoordinateVector::hashCode)
			.sum();
		return result;
	}
	
	@Override
	public boolean isOnFace(final CoordinateVector pos)
	{
		final DistortedCell upStream = getNormalUpstreamCell();
		if (upStream != null)
		{
			return referenceFace.isOnFace(cellsWithTransformation
				                              .get(upStream)
				                              .applyInverse(getNormalUpstreamCell()
					                                            .transformToReferenceCell(
						                                            pos)));
		} else if (getNormalDownstreamCell() != null)
			return referenceFace.isOnFace(cellsWithTransformation
				                              .get(getNormalDownstreamCell())
				                              .applyInverse(getNormalDownstreamCell()
					                                            .transformToReferenceCell(
						                                            pos)));
		else throw new IllegalStateException("Face is isolated");
	}
	
	@Override
	public List<DistortedFace> refine(final Multimap<DistortedCell, DistortedCell> cellMap)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int compareTo(@NotNull final DistortedFace o)
	{
		
		if (o.getDimension() < getDimension())
			return -1;
		if (o.getDimension() > getDimension())
			return 1;
		if (this.equals(o))
			return 0;
		for (int i = 0; i < vertices.length; i++)
			if (vertices[i].compareTo(o.vertices[i]) != 0)
				return vertices[i].compareTo(o.vertices[i]);
		return 0;
	}
	
	@Override
	public DistortedFace getReferenceFace()
	{
		if (getDimension() == 2)
		{
			final CoordinateVector[] vertices = new CoordinateVector[2];
			vertices[0] = CoordinateVector.fromValues(0, 0);
			vertices[1] = CoordinateVector.fromValues(1, 0);
			return new DistortedFace(vertices, isBoundaryFace);
		}
		if (getDimension() == 3)
		{
			final CoordinateVector[] vertices = new CoordinateVector[4];
			vertices[0] = CoordinateVector.fromValues(0, 0, 0);
			vertices[1] = CoordinateVector.fromValues(1, 0, 0);
			vertices[2] = CoordinateVector.fromValues(1, 1, 0);
			vertices[3] = CoordinateVector.fromValues(0, 1, 0);
			return new DistortedFace(vertices, isBoundaryFace);
		}
		throw new IllegalStateException("Dimension must be less than 3");
	}
}
