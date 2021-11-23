package linalg;

import basic.PerformanceArguments;

public class CoordinateDenseMatrix
	extends DenseMatrix
	implements CoordinateMatrix
{
	public CoordinateDenseMatrix(final int i, final int j)
	{
		super(i, j);
//		if (PerformanceArguments.getInstance().executeChecks)
//			if (i > 3 || j > 3)
//				throw new IllegalArgumentException("only 1D, 2D and 3D supported");
	}
	
	public CoordinateDenseMatrix(final Matrix matrix)
	{
		super(matrix);
		if (PerformanceArguments.getInstance().executeChecks)
			if (matrix.getShape()
			          .get(0) > 3 || matrix.getShape()
			                               .get(1) > 3)
				throw new IllegalArgumentException("only 1D, 2D and 3D supported");
	}
	
	public CoordinateDenseMatrix(final double[][] matrix)
	{
		super(matrix);
	}
	
	public static CoordinateDenseMatrix fromValues(final int rows, final int cols, final double... vals)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (rows * cols != vals.length) throw new IllegalArgumentException("matrix does not fit size");
		final CoordinateDenseMatrix ret = new CoordinateDenseMatrix(rows, cols);
		for (int i = 0; i < rows * cols; i++)
		{
			ret.entries[i / cols][i % cols] = vals[i];
		}
		return ret;
	}
	
	public static CoordinateDenseMatrix identity(final int n)
	{
		final CoordinateDenseMatrix ret = new CoordinateDenseMatrix(n, n);
		for (int i = 0; i < n; i++)
			ret.add(1, i, i);
		return ret;
	}
	
	@Override
	public CoordinateVector mvMul(final Vector vector)
	{
		if (entries.length == 2 && entries[0].length == 2)
		{
			final double v0 = vector.at(0);
			final double v1 = vector.at(1);
			return CoordinateVector.fromValues(entries[0][0] * v0 + entries[0][1] * v1,
			                                   entries[1][0] * v0 + entries[1][1] * v1);
		}
		return new CoordinateVector(super.mvMul(vector));
	}
	
	@Override
	public CoordinateVector tvMul(final Vector vector)
	{
		if (entries.length == 2 && entries[0].length == 2)
		{
			final double v0 = vector.at(0);
			final double v1 = vector.at(1);
			return CoordinateVector.fromValues(entries[0][0] * v0 + entries[1][0] * v1,
			                                   entries[0][1] * v0 + entries[1][1] * v1);
		}
		return new CoordinateVector(super.tvMul(vector));
	}
	
	@Override
	public CoordinateDenseMatrix mmMul(final Matrix matrix)
	{
		if (entries.length == 2 && entries[0].length == 2 && matrix.getRows() == 2 && matrix.getCols() == 2)
		{
			return CoordinateDenseMatrix.fromValues(2, 2,
			                                        at(0, 0) * matrix.at(0, 0) + at(0, 1) * matrix.at(1, 0),
			                                        at(0, 0) * matrix.at(0, 1) + at(0, 1) * matrix.at(1, 1),
			                                        at(1, 0) * matrix.at(0, 0) + at(1, 1) * matrix.at(1, 0),
			                                        at(1, 0) * matrix.at(0, 1) + at(1, 1) * matrix.at(1,
			                                                                                          1));
		}
		return new CoordinateDenseMatrix(super.mmMul(matrix));
	}
	
	@Override
	public CoordinateDenseMatrix mtMul(final Matrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (matrix.getCols())) throw new IllegalArgumentException("Incompatible sizes");
		return mmMul(matrix.transpose());
	}
	
	@Override
	public CoordinateDenseMatrix tmMul(final Matrix matrix)
	{
		return new CoordinateDenseMatrix(super.tmMul(matrix));
	}
	
	@Override
	public CoordinateDenseMatrix add(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		final CoordinateDenseMatrix ret = new CoordinateDenseMatrix(this.getCols(), this.getRows());
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.add(other.at(i, j) + at(i, j), i, j);
		return ret;
	}
	
	@Override
	public CoordinateDenseMatrix sub(final Tensor other)
	{
		
		if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		final CoordinateDenseMatrix ret = new CoordinateDenseMatrix(this);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.add(-other.at(i, j), i, j);
		return ret;
	}
	
	@Override
	public CoordinateDenseMatrix inverse()
	{
		if (entries.length == 2)
		{
			final double det = entries[0][0] * entries[1][1] - entries[1][0] * entries[0][1];
			final CoordinateDenseMatrix ret = new CoordinateDenseMatrix(2, 2);
			ret.entries[0][0] = entries[1][1] / det;
			ret.entries[0][1] = -entries[0][1] / det;
			ret.entries[1][0] = -entries[1][0] / det;
			ret.entries[1][1] = entries[0][0] / det;
			return ret;
		}
		return new CoordinateDenseMatrix(super.inverse());
	}
	
	@Override
	public CoordinateDenseMatrix mul(final double scalar)
	{
		final CoordinateDenseMatrix ret = new CoordinateDenseMatrix(getCols(), getRows());
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(scalar * at(i, j), i, j);
		return ret;
	}
	
	@Override
	public CoordinateDenseMatrix transpose()
	{
		final CoordinateDenseMatrix ret = new CoordinateDenseMatrix(entries[0].length, entries.length);
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
	
	@Override
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
