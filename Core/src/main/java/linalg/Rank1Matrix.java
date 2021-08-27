package linalg;

import java.util.List;

public class Rank1Matrix implements Matrix
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
		return List.of(ver, hor);
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return hor.getLength() * ver.getLength();
	}
	
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
	public DenseMatrix mmMul(final Matrix matrix)
	{
		return new DenseMatrix(this).mmMul(matrix);
	}
}
