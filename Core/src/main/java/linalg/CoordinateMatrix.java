package linalg;

public interface CoordinateMatrix extends Matrix
{
	@Override
	CoordinateVector mvMul(final Vector vector);
	
	@Override
	CoordinateVector tvMul(final Vector vector);
	
	@Override
	CoordinateMatrix mmMul(final Matrix matrix);
	
	@Override
	CoordinateMatrix mtMul(final Matrix matrix);
	
	@Override
	public CoordinateMatrix tmMul(final Matrix matrix);
	
	@Override
	public CoordinateMatrix add(final Tensor other);
	
	@Override
	public CoordinateMatrix sub(final Tensor other);
	
	@Override
	public CoordinateMatrix mul(final double scalar);
	
	@Override
	public CoordinateMatrix transpose();
	
	default double determinant()
	{
		if (this.getCols() == 1) return at(0, 0);
		if (this.getCols() == 2) return at(0, 0) * at(1, 1) - at(1, 0) * at(0, 1);
		if (this.getCols() == 3)
			return at(0, 0) * at(1, 1) * at(2, 2) + at(0, 1) * at(1, 2) * at(2, 0) + at(0, 2) * at(1,
			                                                                                       0) * at(
				2, 1) - at(0, 0) * at(1, 2) * at(2, 1) - at(1, 0) * at(2, 2) * at(0, 1) - at(2, 0) * at(
				0, 2) * at(1, 1);
		throw new IllegalStateException("Dimension not allowed");
	}
}
