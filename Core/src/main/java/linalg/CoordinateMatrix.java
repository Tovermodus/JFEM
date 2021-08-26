package linalg;

import basic.PerformanceArguments;

public class CoordinateMatrix extends DenseMatrix
{
	public CoordinateMatrix(final int i, final int j)
	{
		super(i, j);
//		if (PerformanceArguments.getInstance().executeChecks)
//			if (i > 3 || j > 3)
//				throw new IllegalArgumentException("only 1D, 2D and 3D supported");
	}
	
	public CoordinateMatrix(final Matrix matrix)
	{
		super(matrix);
		if (PerformanceArguments.getInstance().executeChecks)
			if (matrix.getShape().get(0) > 3 || matrix.getShape().get(1) > 3)
				throw new IllegalArgumentException("only 1D, 2D and 3D supported");
	}
	
	public CoordinateMatrix(final double[][] matrix)
	{
		super(matrix);
	}
	
	public static CoordinateMatrix fromValues(final int rows, final int cols, final double... vals)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (rows * cols != vals.length) throw new IllegalArgumentException("matrix does not fit size");
		final CoordinateMatrix ret = new CoordinateMatrix(rows, cols);
		for (int i = 0; i < rows * cols; i++)
		{
			ret.entries[i / cols][i % cols] = vals[i];
		}
		return ret;
	}
	
	@Override
	public CoordinateVector mvMul(final Vector vector)
	{
		return new CoordinateVector(super.mvMul(vector));
	}
	
	@Override
	public CoordinateVector tvMul(final Vector vector)
	{
		return new CoordinateVector(super.tvMul(vector));
	}
	
	@Override
	public CoordinateMatrix mmMul(final Matrix matrix)
	{
		return new CoordinateMatrix(super.mmMul(matrix));
	}
	
	@Override
	public CoordinateMatrix mtMul(Matrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (matrix.getCols()))
				throw new IllegalArgumentException("Incompatible sizes");
		return mmMul(matrix.transpose());
	}
	@Override
	public CoordinateMatrix tmMul(final Matrix matrix)
	{
		return new CoordinateMatrix(super.tmMul(matrix));
	}
	
	@Override
	public CoordinateMatrix add(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		final CoordinateMatrix ret = new CoordinateMatrix(this.getCols(), this.getRows());
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.add(other.at(i, j) + at(i, j), i, j);
		return ret;
	}
	
	@Override
	public CoordinateMatrix sub(final Tensor other)
	{
		
		if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		final CoordinateMatrix ret = new CoordinateMatrix(this);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.add(-other.at(i, j), i, j);
		return ret;
	}
	
	@Override
	public CoordinateMatrix inverse()
	{
		if (entries.length == 2)
		{
			final double det = entries[0][0] * entries[1][1] - entries[1][0] * entries[0][1];
			final CoordinateMatrix ret = new CoordinateMatrix(2, 2);
			ret.entries[0][0] = entries[1][1] / det;
			ret.entries[0][1] = -entries[0][1] / det;
			ret.entries[1][0] = -entries[1][0] / det;
			ret.entries[1][1] = entries[0][0] / det;
			return ret;
		}
		return new CoordinateMatrix(super.inverse());
	}
	
	@Override
	public CoordinateMatrix mul(final double scalar)
	{
		final CoordinateMatrix ret = new CoordinateMatrix(getCols(), getRows());
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(scalar * at(i, j), i, j);
		return ret;
	}
	
	@Override
	public CoordinateMatrix transpose()
	{
		final CoordinateMatrix ret = new CoordinateMatrix(entries[0].length, entries.length);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(at(j, i), i, j);
		return ret;
	}
	
	@Override
	public CoordinateVector solve(final Vector rhs)
	{
		if (entries.length == 2) return inverse().mvMul(rhs);
		return new CoordinateVector(super.solve(rhs));
	}
	
	public double determinant()
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
