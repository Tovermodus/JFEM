package distorted.geometry;
import basic.Cell;
import basic.CellWithReferenceCell;
import basic.VectorFunction;
import com.google.common.collect.ImmutableSet;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Newton;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;

import java.util.Arrays;
import java.util.List;

public class DistortedCell implements CellWithReferenceCell<DistortedCell, DistortedFace>
{
	TPCell referenceCell;
	CoordinateVector[] transformationCoefficients; //in order xyz, xy, xz, yz, x, y, z, 1
	public DistortedCell(CoordinateVector[] vertices)
	{
		transformationCoefficients = new CoordinateVector[4];
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
	@Override
	public int getDimension()
	{
		return 2;
	}
	
	@Override
	public ImmutableSet<DistortedFace> getFaces()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public boolean isInCell(CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public CoordinateVector center()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public VectorFunction getOuterNormal(DistortedFace face)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public List<DistortedCell> refine(List<DistortedFace> refinedFaces)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int compareTo(@NotNull DistortedCell o)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public DistortedCell getReferenceCell()
	{
		throw new UnsupportedOperationException("not implemented yet");
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
		return Newton.solve(CoordinateVector.repeat(0.5, getDimension()), pos,
			getTransformationFromReferenceCell());
	}
	
	@Override
	public CoordinateVector transformFromReferenceCell(CoordinateVector pos)
	{
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
				.add(transformationCoefficients[3].mul(pos.x()*pos.z()))
				.add(transformationCoefficients[4].mul(pos.y()*pos.z()))
				.add(transformationCoefficients[5].mul(pos.x()))
				.add(transformationCoefficients[6].mul(pos.y()))
				.add(transformationCoefficients[7].mul(pos.z()))
				.add(transformationCoefficients[8]);
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
}
