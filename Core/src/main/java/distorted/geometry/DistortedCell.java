package distorted.geometry;
import basic.CellWithReferenceCell;
import basic.PerformanceArguments;
import basic.VectorFunction;
import com.google.common.collect.ImmutableSet;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Newton;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class DistortedCell implements CellWithReferenceCell<DistortedCell, DistortedFace>
{
	
	final TPCell referenceCell;
	private final CoordinateVector[] transformationCoefficients; //in order xyz, xy, xz, yz, x, y, z, 1
	private final int dimension;
	Set<DistortedFace> faces;
	public DistortedCell(CoordinateVector[] vertices)
	{
		dimension = vertices[0].getLength();
		referenceCell = TPCell.unitHyperCube(dimension);
		faces = new HashSet<>(2*dimension);
		if(PerformanceArguments.getInstance().executeChecks)
		{
			if (vertices.length != (int)Math.pow(2,dimension))
				throw new IllegalArgumentException("Wrong number of vertices");
		}
		transformationCoefficients = new CoordinateVector[(int)Math.pow(2,dimension)];
		if(dimension == 1)
		{
			transformationCoefficients[0] = vertices[1]
				.sub(vertices[0]);
			transformationCoefficients[1] = vertices[0];
		}
		if(dimension == 2)
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
		if(dimension == 3)
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
	}
	@Override
	public int getDimension()
	{
		return dimension;
	}
	
	@Override
	public ImmutableSet<DistortedFace> getFaces()
	{
		return ImmutableSet.copyOf(faces);
	}
	
	@Override
	public boolean isInCell(CoordinateVector pos)
	{
		return referenceCell.isInCell(transformToReferenceCell(pos));
	}
	
	@Override
	public CoordinateVector center()
	{
		return transformFromReferenceCell(referenceCell.center());
	}
	
	@Override
	public VectorFunction getOuterNormal(DistortedFace face)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!faces.contains(face))
				throw new IllegalArgumentException("face does not belong to cell");
		boolean invertNormal = (center().sub(face.center())).inner(face.getNormal().value(face.center())) > 0;
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
				public CoordinateVector value(CoordinateVector pos)
				{
					return face.getNormal().value(pos).mul(-1);
				}
			};
		else return face.getNormal();
	}
	
	@Override
	public List<DistortedCell> refine(List<DistortedFace> refinedFaces)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int compareTo(@NotNull DistortedCell o)
	{
		if (o.getDimension() < getDimension())
			return -1;
		if (o.getDimension() > getDimension())
			return 1;
		for(int i = 0; i < transformationCoefficients.length; i++)
			if(transformationCoefficients[i].compareTo(o.transformationCoefficients[i]) != 0)
				return transformationCoefficients[i].compareTo(o.transformationCoefficients[i]);
		return 0;
	}
	
	@Override
	public DistortedCell getReferenceCell()
	{
		if(getDimension() == 1)
		{
			CoordinateVector [] vertices = new CoordinateVector[2];
			vertices[0] = CoordinateVector.fromValues(0);
			vertices[1] = CoordinateVector.fromValues(1);
			return new DistortedCell(vertices);
		}
		if(getDimension() == 2)
		{
			CoordinateVector [] vertices = new CoordinateVector[4];
			vertices[0] = CoordinateVector.fromValues(0,0);
			vertices[1] = CoordinateVector.fromValues(1,0);
			vertices[2] = CoordinateVector.fromValues(1,1);
			vertices[3] = CoordinateVector.fromValues(0,1);
			return new DistortedCell(vertices);
		}
		if(getDimension() == 3)
		{
			CoordinateVector [] vertices = new CoordinateVector[8];
			vertices[0] = CoordinateVector.fromValues(0,0,0);
			vertices[1] = CoordinateVector.fromValues(1,0,0);
			vertices[2] = CoordinateVector.fromValues(1,1,0);
			vertices[3] = CoordinateVector.fromValues(0,1,0);
			vertices[4] = CoordinateVector.fromValues(0,0,1);
			vertices[5] = CoordinateVector.fromValues(1,0,1);
			vertices[6] = CoordinateVector.fromValues(1,1,1);
			vertices[7] = CoordinateVector.fromValues(0,1,1);
			return new DistortedCell(vertices);
		}
		throw new IllegalStateException("Dimension must be less than 3");
	}
	
	@Override
	public CoordinateMatrix transformationGradientFromReferenceCell(CoordinateVector pos)
	{
		if(getDimension() == 1)
		{
			return CoordinateMatrix.fromValues(1,1,
				transformationCoefficients[0].at(0));
		}
		if(getDimension() == 2)
		{
			CoordinateMatrix ret = new CoordinateMatrix(2,2);
			ret.addColumn(transformationCoefficients[0].mul(pos.y())
				.add(transformationCoefficients[1]), 0);
			ret.addColumn(transformationCoefficients[0].mul(pos.x())
				.add(transformationCoefficients[2]), 1);
			return ret;
		}
		if(getDimension() == 3)
		{
			CoordinateMatrix ret = new CoordinateMatrix(3,3);
			ret.addColumn(transformationCoefficients[0].mul(pos.y()*pos.z())
				.add(transformationCoefficients[1].mul(pos.y()))
				.add(transformationCoefficients[2].mul(pos.z()))
				.add(transformationCoefficients[4]), 0);
			ret.addColumn(transformationCoefficients[0].mul(pos.x()*pos.z())
				.add(transformationCoefficients[1].mul(pos.x()))
				.add(transformationCoefficients[3].mul(pos.z()))
				.add(transformationCoefficients[5]), 1);
			ret.addColumn(transformationCoefficients[0].mul(pos.y()*pos.x())
				.add(transformationCoefficients[2].mul(pos.x()))
				.add(transformationCoefficients[3].mul(pos.y()))
				.add(transformationCoefficients[6]), 2);
			return ret;
		}
		throw new IllegalStateException("Dimension must be less than 3");
	}
	
	@Override
	public CoordinateMatrix transformationGradientToReferenceCell(CoordinateVector pos)
	{
		return transformationGradientFromReferenceCell(transformToReferenceCell(pos)).inverse();
	}
	
	@Override
	public CoordinateVector transformToReferenceCell(CoordinateVector pos)
	{
		if(PerformanceArguments.getInstance().executeChecks)
		{
			if (pos.getLength() != getDimension())
				throw new IllegalArgumentException("coordinate has wrong dimension");
		}
		return Newton.solve(CoordinateVector.repeat(0.5, getDimension()), pos,
			getTransformationFromReferenceCell());
	}
	
	@Override
	public CoordinateVector transformFromReferenceCell(CoordinateVector pos)
	{
		if(PerformanceArguments.getInstance().executeChecks)
		{
			if (pos.getLength() != getDimension())
				throw new IllegalArgumentException("coordinate has wrong dimension");
		}
		if(getDimension() == 1)
		{
			return transformationCoefficients[0].mul(pos.x()).add(transformationCoefficients[1]);
		}
		if(getDimension() == 2)
		{
			return transformationCoefficients[0].mul(pos.x()*pos.y())
				.add(transformationCoefficients[1].mul(pos.x()))
				.add(transformationCoefficients[2].mul(pos.y()))
				.add(transformationCoefficients[3]);
		}
		if(getDimension() == 3)
		{
			return transformationCoefficients[0].mul(pos.x()*pos.y()*pos.z())
				.add(transformationCoefficients[1].mul(pos.x()*pos.y()))
				.add(transformationCoefficients[2].mul(pos.x()*pos.z()))
				.add(transformationCoefficients[3].mul(pos.y()*pos.z()))
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
			public CoordinateVector value(CoordinateVector pos)
			{
				return transformFromReferenceCell(pos);
			}
			
			@Override
			public CoordinateMatrix gradient(CoordinateVector pos)
			{
				return transformationGradientFromReferenceCell(pos);
			}
		};
	}
	@Override
	public boolean equals(Object o)
	{
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		DistortedCell that = (DistortedCell) o;
		return Arrays.deepEquals(transformationCoefficients, that.transformationCoefficients);
	}
	
	@Override
	public int hashCode()
	{
		return Arrays.deepHashCode(transformationCoefficients);
	}
}
