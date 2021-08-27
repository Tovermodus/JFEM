package linalg;

public class Rank1CoordinateMatrix extends Rank1Matrix implements CoordinateMatrix
{
	@Override
	public CoordinateVector mvMul(final Vector vector)
	{
		return (CoordinateVector) hor.mul(ver.inner(vector));
	}
	
	@Override
	public CoordinateDenseMatrix mmMul(Matrix matrix)
	{
		return new CoordinateDenseMatrix(super.mmMul(matrix));
	}
	
	public Rank1CoordinateMatrix(final CoordinateVector ver, final CoordinateVector hor)
	{
		super(ver, hor);
	}
	
	@Override
	public CoordinateVector tvMul(final Vector vector)
	{
		return (CoordinateVector) super.tvMul(vector);
	}
	
	@Override
	public CoordinateDenseMatrix mtMul(final Matrix matrix)
	{
		return new CoordinateDenseMatrix(super.mtMul(matrix));
	}
	
	@Override
	public CoordinateMatrix tmMul(final Matrix matrix)
	{
		return new CoordinateDenseMatrix(super.tmMul(matrix));
	}
	
	@Override
	public CoordinateMatrix sub(final Tensor other)
	{
		return new CoordinateDenseMatrix(super.sub(other));
	}
	
	@Override
	public CoordinateDenseMatrix add(final Tensor other)
	{
		return new CoordinateDenseMatrix(super.add(other));
	}
	
	@Override
	public Rank1CoordinateMatrix mul(final double scalar)
	{
		return new Rank1CoordinateMatrix((CoordinateVector) ver.mul(scalar), (CoordinateVector) hor);
	}
	
	@Override
	public Rank1CoordinateMatrix transpose()
	{
		return new Rank1CoordinateMatrix((CoordinateVector) hor, (CoordinateVector) ver);
	}
	
	public void deleteRow(final int row)
	{
		ver.set(0, row);
	}
	
	public void deleteColumn(final int column)
	{
		hor.set(0, column);
	}
	
	public void mulInPlace(final double scalar)
	{
		ver.mulInPlace(scalar);
	}
}
