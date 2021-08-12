package distorted.geometry;

import basic.CellWithReferenceCell;
import basic.DoubleCompare;
import basic.PerformanceArguments;
import basic.VectorFunction;
import com.google.common.collect.ImmutableSet;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.Newton;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;

import java.util.*;

public class DistortedCell implements CellWithReferenceCell<DistortedCell, DistortedFace>
{
	final double MAXIMUM_WARP = 1e-10; //25°
	final TPCell referenceCell;
	final CoordinateVector[] vertices;
	private final CoordinateVector[] transformationCoefficients;//in order xyz, xy, xz, yz, x, y, z, 1
	private final int dimension;
	Set<DistortedFace> faces;
	
	public DistortedCell(final CoordinateVector... vertices)
	{
		this.vertices = vertices;
		dimension = vertices[0].getLength();
		referenceCell = TPCell.unitHyperCube(dimension);
		faces = new HashSet<>(2 * dimension);
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (vertices.length != (int) Math.pow(2, dimension))
				throw new IllegalArgumentException("Wrong number of vertices");
		}
		transformationCoefficients = new CoordinateVector[(int) Math.pow(2, dimension)];
		if (getDimension() == 1)
		{
			transformationCoefficients[0] = vertices[1]
				.sub(vertices[0]);
			transformationCoefficients[1] = vertices[0];
		}
		if (getDimension() == 2)
		{
			transformationCoefficients[0] = vertices[0]
				.add(vertices[2])
				.sub(vertices[1])
				.sub(vertices[3]);
			transformationCoefficients[1] = vertices[1]
				.sub(vertices[0]);
			transformationCoefficients[2] = vertices[3]
				.sub(vertices[0]);
			transformationCoefficients[3] = vertices[0];
		}
		if (getDimension() == 3)
		{
			transformationCoefficients[0] = vertices[1]
				.add(vertices[3])
				.add(vertices[4])
				.add(vertices[6])
				.sub(vertices[0])
				.sub(vertices[2])
				.sub(vertices[5])
				.sub(vertices[7]);
			transformationCoefficients[1] = vertices[0]
				.add(vertices[2])
				.sub(vertices[1])
				.sub(vertices[3]);
			transformationCoefficients[2] = vertices[0]
				.add(vertices[5])
				.sub(vertices[4])
				.sub(vertices[1]);
			transformationCoefficients[3] = vertices[0]
				.add(vertices[7])
				.sub(vertices[3])
				.sub(vertices[4]);
			transformationCoefficients[4] = vertices[1]
				.sub(vertices[0]);
			transformationCoefficients[5] = vertices[3]
				.sub(vertices[0]);
			transformationCoefficients[6] = vertices[4]
				.sub(vertices[0]);
			transformationCoefficients[7] = vertices[0];
		}
		if (!verticesHaveCorrectPosition())
			throw new IllegalArgumentException("vertices are in wrong order or hexahedron is warped too " +
				                                   "much");
	}
	
	@Override
	public int getDimension()
	{
		return dimension;
	}
	
	public int getPositionOfVertex(final CoordinateVector vertex)
	{
		for (int i = 0; i < vertices.length; i++)
			if (vertices[i].almostEqualMute(vertex))
				return i;
		return -1;
	}
	
	public DistortedFace getFaceFromVertexNumbers(final int... numbers)
	{
		for (final DistortedFace face : faces)
		{
			System.out.println(face);
			System.out.println(vertices[numbers[0]] + " " + vertices[numbers[1]]);
			if (Arrays.stream(numbers).mapToObj(i -> vertices[i]).allMatch(face::isOnFace))
				return face;
		}
		throw new IllegalArgumentException("Cell has no Face with the numbers" + Arrays.toString(numbers));
	}
	
	public IntCoordinates mapPositionToReferencePosition(final int position)
	{
		if (getDimension() == 1)
		{
			return new IntCoordinates(position);
		}
		if (getDimension() == 2)
		{
			if (position == 0)
				return new IntCoordinates(0, 0);
			if (position == 1)
				return new IntCoordinates(1, 0);
			if (position == 2)
				return new IntCoordinates(1, 1);
			if (position == 3)
				return new IntCoordinates(0, 1);
		}
		if (getDimension() == 3)
		{
			if (position == 0)
				return new IntCoordinates(0, 0, 0);
			if (position == 1)
				return new IntCoordinates(1, 0, 0);
			if (position == 2)
				return new IntCoordinates(1, 1, 0);
			if (position == 3)
				return new IntCoordinates(0, 1, 0);
			if (position == 4)
				return new IntCoordinates(0, 0, 1);
			if (position == 5)
				return new IntCoordinates(1, 0, 1);
			if (position == 6)
				return new IntCoordinates(1, 1, 1);
			if (position == 7)
				return new IntCoordinates(0, 1, 1);
		}
		throw new IllegalArgumentException("position or dimension is wrong");
	}
	
	private boolean verticesHaveCorrectPosition()
	{
		
		final CoordinateVector rating;
		if (dimension == 2)
			rating = CoordinateVector.fromValues(1, 1);
		else
			rating = CoordinateVector.fromValues(1, 1, 2);
		final OptionalDouble minSum = Arrays.stream(vertices).mapToDouble(rating::inner).min();
		double min = 0;
		if (minSum.isPresent())
			min = minSum.getAsDouble();
		else
			throw new IllegalStateException("No vertices");
		if (!DoubleCompare.almostEqual(vertices[0].inner(rating), min))
		{
			System.out.println("min wrong" + min + " " + vertices[0].inner(rating));
			return false;
		}
		if (getDimension() == 2)
			if (!checkQuadrilateral(vertices[0], vertices[3], vertices[2], vertices[1]))
			{
				System.out.println("order wrong");
				return false;
			}
		if (getDimension() == 3)
		{
			if (!checkQuadrilateral(vertices[0], vertices[3], vertices[2], vertices[1]))
			{
				System.out.println("bottom wrong");
				return false;
			}
			if (!checkQuadrilateral(vertices[0], vertices[1], vertices[5], vertices[4]))
			{
				System.out.println("front wrong");
				return false;
			}
			if (!checkQuadrilateral(vertices[0], vertices[4], vertices[7], vertices[3]))
			{
				System.out.println("left wrong");
				return false;
			}
			if (!checkQuadrilateral(vertices[6], vertices[7], vertices[4], vertices[5]))
			{
				System.out.println("top wrong");
				return false;
			}
			if (!checkQuadrilateral(vertices[6], vertices[2], vertices[3], vertices[7]))
			{
				System.out.println("back wrong");
				return false;
			}
			if (!checkQuadrilateral(vertices[6], vertices[5], vertices[1], vertices[2]))
			{
				System.out.println("right wrong");
				return false;
			}
		}
		return true;
	}
	
	private boolean checkQuadrilateral(CoordinateVector c1, CoordinateVector c2, CoordinateVector c3,
	                                   CoordinateVector c4)
	{
		CoordinateVector center = center();
		if (getDimension() == 2)
		{
			c1 = c1.addCoordinate(0);
			c2 = c2.addCoordinate(0);
			c3 = c3.addCoordinate(0);
			c4 = c4.addCoordinate(0);
			center = center.addCoordinate(1);
		}
		final CoordinateVector v1 = c2.sub(c1);
		final CoordinateVector v2 = c3.sub(c2);
		final CoordinateVector v3 = c4.sub(c3);
		final CoordinateVector v4 = c1.sub(c4);
		final CoordinateVector c12 = v1.cross(v2);
		final CoordinateVector c23 = v2.cross(v3);
		final CoordinateVector c34 = v3.cross(v4);
		final CoordinateVector c41 = v4.cross(v1);
		c12.mulInPlace(1. / c12.euclidianNorm());
		c23.mulInPlace(1. / c23.euclidianNorm());
		c34.mulInPlace(1. / c34.euclidianNorm());
		c41.mulInPlace(1. / c41.euclidianNorm());
		final CoordinateVector mean = c1.add(c2).add(c3).add(c4).mul(1. / 4);
		if (c12.inner(mean.sub(center)) <= 0)
			return false;
		if (c12.inner(c23) < 1 - MAXIMUM_WARP)
			return false;
		if (c12.inner(c34) < 1 - MAXIMUM_WARP)
			return false;
		if (c12.inner(c41) < 1 - MAXIMUM_WARP)
			return false;
		if (c23.inner(c34) < 1 - MAXIMUM_WARP)
			return false;
		if (c23.inner(c41) < 1 - MAXIMUM_WARP)
			return false;
		if (c34.inner(c41) < 1 - MAXIMUM_WARP)
			return false;
		return true;
	}
	
	@Override
	public ImmutableSet<DistortedFace> getFaces()
	{
		return ImmutableSet.copyOf(faces);
	}
	
	@Override
	public boolean isInCell(final CoordinateVector pos)
	{
		return referenceCell.isInCell(transformToReferenceCell(pos));
	}
	
	@Override
	public CoordinateVector center()
	{
		return Arrays.stream(vertices).reduce(new CoordinateVector(dimension), CoordinateVector::add).mul(
			1. / vertices.length);//transformFromReferenceCell(referenceCell.center());
	}
	
	@Override
	public VectorFunction getOuterNormal(final DistortedFace face)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!faces.contains(face))
				throw new IllegalArgumentException("face does not belong to cell");
		final boolean invertNormal = (center().sub(face.center())).inner(
			face.getNormal().value(face.center())) > 0;
		if (invertNormal)
			return new VectorFunction()
			{
				@Override
				public int getRangeDimension()
				{
					return face.getNormal().getRangeDimension();
				}
				
				@Override
				public int getDomainDimension()
				{
					return face.getNormal().getDomainDimension();
				}
				
				@Override
				public CoordinateVector value(final CoordinateVector pos)
				{
					return face.getNormal().value(pos).mul(-1);
				}
			};
		else return face.getNormal();
	}
	
	@Override
	public List<DistortedCell> refine(final List<DistortedFace> refinedFaces)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int compareTo(@NotNull final DistortedCell o)
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
	public DistortedCell getReferenceCell()
	{
		if (getDimension() == 1)
		{
			final CoordinateVector[] vertices = new CoordinateVector[2];
			vertices[0] = CoordinateVector.fromValues(0);
			vertices[1] = CoordinateVector.fromValues(1);
			return new DistortedCell(vertices);
		}
		if (getDimension() == 2)
		{
			final CoordinateVector[] vertices = new CoordinateVector[4];
			vertices[0] = CoordinateVector.fromValues(0, 0);
			vertices[1] = CoordinateVector.fromValues(1, 0);
			vertices[2] = CoordinateVector.fromValues(1, 1);
			vertices[3] = CoordinateVector.fromValues(0, 1);
			return new DistortedCell(vertices);
		}
		if (getDimension() == 3)
		{
			final CoordinateVector[] vertices = new CoordinateVector[8];
			vertices[0] = CoordinateVector.fromValues(0, 0, 0);
			vertices[1] = CoordinateVector.fromValues(1, 0, 0);
			vertices[2] = CoordinateVector.fromValues(1, 1, 0);
			vertices[3] = CoordinateVector.fromValues(0, 1, 0);
			vertices[4] = CoordinateVector.fromValues(0, 0, 1);
			vertices[5] = CoordinateVector.fromValues(1, 0, 1);
			vertices[6] = CoordinateVector.fromValues(1, 1, 1);
			vertices[7] = CoordinateVector.fromValues(0, 1, 1);
			return new DistortedCell(vertices);
		}
		throw new IllegalStateException("Dimension must be less than 3");
	}
	
	@Override
	public CoordinateMatrix transformationGradientFromReferenceCell(final CoordinateVector pos)
	{
		if (getDimension() == 1)
		{
			return CoordinateMatrix.fromValues(1, 1,
			                                   transformationCoefficients[0].at(0));
		}
		if (getDimension() == 2)
		{
			final CoordinateMatrix ret = new CoordinateMatrix(2, 2);
			ret.addColumn(transformationCoefficients[0].mul(pos.y())
			                                           .add(transformationCoefficients[1]), 0);
			ret.addColumn(transformationCoefficients[0].mul(pos.x())
			                                           .add(transformationCoefficients[2]), 1);
			return ret;
		}
		if (getDimension() == 3)
		{
			final CoordinateMatrix ret = new CoordinateMatrix(3, 3);
			ret.addColumn(transformationCoefficients[0].mul(pos.y() * pos.z())
			                                           .add(transformationCoefficients[1].mul(pos.y()))
			                                           .add(transformationCoefficients[2].mul(pos.z()))
			                                           .add(transformationCoefficients[4]), 0);
			ret.addColumn(transformationCoefficients[0].mul(pos.x() * pos.z())
			                                           .add(transformationCoefficients[1].mul(pos.x()))
			                                           .add(transformationCoefficients[3].mul(pos.z()))
			                                           .add(transformationCoefficients[5]), 1);
			ret.addColumn(transformationCoefficients[0].mul(pos.y() * pos.x())
			                                           .add(transformationCoefficients[2].mul(pos.x()))
			                                           .add(transformationCoefficients[3].mul(pos.y()))
			                                           .add(transformationCoefficients[6]), 2);
			return ret;
		}
		throw new IllegalStateException("Dimension must be less than 3");
	}
	
	@Override
	public CoordinateMatrix transformationGradientToReferenceCell(final CoordinateVector pos)
	{
		return transformationGradientFromReferenceCell(transformToReferenceCell(pos)).inverse();
	}
	
	@Override
	public CoordinateVector transformToReferenceCell(final CoordinateVector pos)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (pos.getLength() != getDimension())
				throw new IllegalArgumentException("coordinate has wrong dimension");
		}
		return Newton.solve(CoordinateVector.repeat(0.5, getDimension()), pos,
		                    getTransformationFromReferenceCell());
	}
	
	@Override
	public CoordinateVector transformFromReferenceCell(final CoordinateVector pos)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (pos.getLength() != getDimension())
				throw new IllegalArgumentException("coordinate has wrong dimension");
		}
		if (getDimension() == 1)
		{
			return transformationCoefficients[0].mul(pos.x()).add(transformationCoefficients[1]);
		}
		if (getDimension() == 2)
		{
			return transformationCoefficients[0].mul(pos.x() * pos.y())
			                                    .add(transformationCoefficients[1].mul(pos.x()))
			                                    .add(transformationCoefficients[2].mul(pos.y()))
			                                    .add(transformationCoefficients[3]);
		}
		if (getDimension() == 3)
		{
			return transformationCoefficients[0].mul(pos.x() * pos.y() * pos.z())
			                                    .add(transformationCoefficients[1].mul(pos.x() * pos.y()))
			                                    .add(transformationCoefficients[2].mul(pos.x() * pos.z()))
			                                    .add(transformationCoefficients[3].mul(pos.y() * pos.z()))
			                                    .add(transformationCoefficients[4].mul(pos.x()))
			                                    .add(transformationCoefficients[5].mul(pos.y()))
			                                    .add(transformationCoefficients[6].mul(pos.z()))
			                                    .add(transformationCoefficients[7]);
		}
		throw new IllegalStateException("Dimension must be less than 3");
	}
	
	public VectorFunction getTransformationFromReferenceCell()
	{
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return getDimension();
			}
			
			@Override
			public int getDomainDimension()
			{
				return getDimension();
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return transformFromReferenceCell(pos);
			}
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return transformationGradientFromReferenceCell(pos);
			}
		};
	}
	
	@Override
	public boolean equals(final Object o)
	{
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		final DistortedCell that = (DistortedCell) o;
		for (int i = 0; i < vertices.length; i++)
		{
			boolean otherHasVertex = false;
			for (int j = 0; j < vertices.length; j++)
				if (vertices[i].almostEqualMute(that.vertices[j]))
					otherHasVertex = true;
			if (!otherHasVertex)
				return false;
		}
		return true;
	}
	
	@Override
	public int hashCode()
	{
		return Arrays.stream(vertices).mapToInt(CoordinateVector::hashCode).sum();
	}
	
	public List<CoordinateVector[]> getSides()
	{
		if (getDimension() == 2)
		{
			final List<CoordinateVector[]> sides = new ArrayList<>(4);
			sides.add(new CoordinateVector[]{vertices[0], vertices[1]});
			sides.add(new CoordinateVector[]{vertices[1], vertices[2]});
			sides.add(new CoordinateVector[]{vertices[2], vertices[3]});
			sides.add(new CoordinateVector[]{vertices[3], vertices[0]});
			return sides;
		}
		if (getDimension() == 3)
		{
			final List<CoordinateVector[]> sides = new ArrayList<>(6);
			sides.add(new CoordinateVector[]{vertices[0], vertices[3], vertices[2], vertices[1]});
			sides.add(new CoordinateVector[]{vertices[0], vertices[1], vertices[5], vertices[4]});
			sides.add(new CoordinateVector[]{vertices[0], vertices[4], vertices[7], vertices[3]});
			sides.add(new CoordinateVector[]{vertices[6], vertices[7], vertices[4], vertices[5]});
			sides.add(new CoordinateVector[]{vertices[6], vertices[2], vertices[3], vertices[7]});
			sides.add(new CoordinateVector[]{vertices[6], vertices[5], vertices[1], vertices[2]});
			return sides;
		} else throw new IllegalStateException("Only Dimension 2 and 3 allowed");
	}
	
	@Override
	public String toString()
	{
		return "DistortedCell{" +
			"vertices=" + Arrays.toString(vertices) +
			'}';
	}
}
