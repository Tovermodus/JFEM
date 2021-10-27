package linalg;

import basic.PerformanceArguments;

import java.util.ArrayList;
import java.util.List;

public class Rank1Matrix
	implements Matrix
{
	public final MutableVector hor;
	public final MutableVector ver;
	
	public Rank1Matrix(final MutableVector ver, final MutableVector hor)
	{
		this.ver = ver;
		this.hor = hor;
	}
	
	public Rank1Matrix(final Vector ver, final Vector hor)
	{
		this.ver = new DenseVector(ver);
		this.hor = new DenseVector(hor);
	}
	
	@Override
	public double at(final int... coordinates)
	{
		return ver.at(coordinates[0]) * hor.at(coordinates[1]);
	}
	
	@Override
	public DenseMatrix add(final Tensor other)
	{
		return new DenseMatrix(this).add(other);
	}
	
	@Override
	public Rank1Matrix mul(final double scalar)
	{
		return new Rank1Matrix(ver.mul(scalar), hor);
	}
	
	@Override
	public List<Vector> unfoldDimension(final int dimension)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (dimension > 1) throw new IllegalArgumentException("Matrix is two dimensional");
		final List<Vector> ret = new ArrayList<>();
		if (dimension == 0)
		{
			for (int i = 0; i < getRows(); i++)
			{
				
				ret.add(hor.mul(ver.at(i)));
			}
		}
		if (dimension == 1)
		{
			for (int i = 0; i < getCols(); i++)
			{
				
				ret.add(ver.mul(hor.at(i)));
			}
		}
		return ret;
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return hor.getLength() * ver.getLength();
	}
	
	public DenseMatrix toDense()
	{
		return new DenseMatrix(this);
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (!(obj instanceof Matrix)) return false;
		return almostEqual((Tensor) obj);
	}
//	@Override
//	public double frobeniusInner(final Matrix other)
//	{
//		return ver.inner(other.mvMul(hor));
//	}
	
	@Override
	public boolean isSparse()
	{
		return false;
	}
	
	@Override
	public IntCoordinates getShape()
	{
		return new IntCoordinates(ver.getLength(), hor.getLength());
	}
	
	@Override
	public Rank1Matrix transpose()
	{
		return new Rank1Matrix(hor, ver);
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		return ver.mul(hor.inner(vector));
	}
	
	@Override
	public Rank1Matrix mmMul(final Matrix matrix)
	{
		
		final Vector tvm = matrix.tvMul(hor);
		return new Rank1Matrix(ver, tvm);
	}
}
