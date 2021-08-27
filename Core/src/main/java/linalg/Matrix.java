package linalg;

import basic.PerformanceArguments;

import java.util.List;

public interface Matrix extends Tensor, VectorMultiplyable
{
	@Override
	default int getTVectorSize()
	{
		return getCols();
	}
	
	@Override
	Matrix add(Tensor other);
	
	@Override
	default int getVectorSize()
	{
		return getRows();
	}
	
	@Override
	default Matrix sub(final Tensor other)
	{
		if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Matrices are of different size");
		return add(other.mul(-1.));
	}
	
	@Override
	default SparseMatrix slice(final IntCoordinates start, final IntCoordinates end)
	{
		final SparseMatrix ret = new SparseMatrix(end.get(0) - start.get(0), end.get(1) - start.get(1));
		for (final IntCoordinates c : new IntCoordinates.Range(start, end))
		{
			ret.set(at(c), c.sub(start));
		}
		return ret;
	}
	
	default DenseVector getColumn(final int col)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (col < 0 || col > getCols()) throw new IllegalArgumentException("col out of bounds");
		final DenseVector ret = new DenseVector(getRows());
		for (int i = 0; i < getRows(); i++)
			ret.set(at(i, col), i);
		return ret;
	}
	
	default DenseVector getRow(final int row)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (row < 0 || row > getRows()) throw new IllegalArgumentException("row out of bounds");
		final DenseVector ret = new DenseVector(getCols());
		for (int i = 0; i < getCols(); i++)
			ret.set(at(row, i), i);
		return ret;
	}
	
	@Override
	Matrix mul(double scalar);
	
	default int getRows()
	{
		return getShape().get(0);
	}
	
	default int getCols()
	{
		return getShape().get(1);
	}
	
	@Override
	List<Vector> unfoldDimension(int dimension);
	
	@Override
	Matrix transpose();
	
	@Override
	Vector mvMul(Vector vector);
	
	@Override
	default Vector tvMul(final Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != vector.getLength()) throw new IllegalArgumentException("Incompatible sizes");
		return transpose().mvMul(vector);
	}
	
	default double frobeniusInner(final Matrix other)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible " + "sizes");
		return getShape().range().stream().mapToDouble(c -> at(c) * other.at(c)).sum();
	}
	
	Matrix mmMul(Matrix matrix);
	
	default Matrix tmMul(final Matrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != (matrix.getRows())) throw new IllegalArgumentException("Incompatible sizes");
		return transpose().mmMul(matrix);
	}
	
	default Matrix mtMul(final Matrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (matrix.getCols())) throw new IllegalArgumentException("Incompatible sizes");
		return mmMul(matrix.transpose());
	}
	
	default SparseMatrix diag()
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != getCols()) throw new IllegalArgumentException("Only for aquare matrices");
		final SparseMatrix ret = new SparseMatrix(getRows(), getCols());
		for (int i = 0; i < getRows(); i++)
		{
			ret.add(at(i, i), i, i);
		}
		return ret;
	}
	
	@Override
	default String printFormatted(double... tol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (tol.length > 1) throw new IllegalArgumentException("Only one Tolerance accepted");
		if (tol.length == 0) tol = new double[]{1e-14};
		String ret = "[";
		for (int i = 0; i < getRows(); i++)
		{
			for (int j = 0; j < getCols(); j++)
			{
				if (Math.abs(at(i, j)) < tol[0]) ret = ret.concat("<>........  ");
				else
				{
					if (at(i, j) >= 0) ret = ret.concat("+");
					ret = ret.concat(String.format("%6.3e", at(i, j)) + "  ");
				}
			}
			ret = ret.concat("\n ");
		}
		ret = ret.concat("]");
		return ret;
	}
}
