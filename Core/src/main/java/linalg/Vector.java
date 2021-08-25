package linalg;

import basic.PerformanceArguments;

import java.util.List;

public interface Vector extends Tensor
{
	default int getLength()
	{
		return getShape().size();
	}
	
	@Override
	default List<Tensor> unfoldDimension(final int dimension)
	{
		throw new UnsupportedOperationException("Vector can't unfold");
	}
	
	@Override
	Vector add(Tensor other);
	
	default double manhattanNorm()
	{
		return getShape().range().stream().mapToDouble(this::at).map(Math::abs).sum();
	}
	
	@Override
	default Vector sub(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Vectors are of different size");
		return add(other.mul(-1.));
	}
	
	@Override
	default DenseVector slice(final IntCoordinates start, final IntCoordinates end)
	{
		final DenseVector ret = new DenseVector(end.get(0) - start.get(0));
		int i = 0;
		for (final IntCoordinates c : new IntCoordinates.Range(start, end))
		{
			ret.set(at(c), i++);
		}
		return ret;
	}
	
	default DenseVector slice(final int start, final int end)
	{
		final DenseVector ret = new DenseVector(end - start);
		int i = 0;
		for (int c = start; c < end; c++)
		{
			ret.set(at(c), i++);
		}
		return ret;
	}
	
	default SparseMatrix asMatrix()
	{
		final SparseMatrix ret = new SparseMatrix(getLength(), 1);
		for (int i = 0; i < getLength(); i++)
		{
			ret.add(at(i), i, 0);
		}
		return ret;
	}
	
	@Override
	Vector mul(double scalar);
	
	Matrix outer(Vector other);
	
	default double inner(final Vector other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getLength() != other.getLength())
				throw new IllegalArgumentException("Vectors are of different size");
		return getShape().range().stream().mapToDouble(c -> at(c) * other.at(c)).sum();
	}
	
	default double euclidianNorm()
	{
		return Math.sqrt(this.inner(this));
	}
	
	default Vector vmMul(final Matrix matrix)
	{
		return matrix.tvMul(this);
	}
	
	@Override
	default String printFormatted(double... tol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (tol.length > 1)
				throw new IllegalArgumentException("Only one Tolerance accepted");
		if (tol.length == 0)
			tol = new double[]{1e-14};
		String ret = "[";
		for (int i = 0; i < getLength(); i++)
		{
			if (Math.abs(at(i)) < tol[0])
				ret = ret.concat("<>........  ");
			else
			{
				if (at(i) >= 0)
					ret = ret.concat("+");
				ret = ret.concat(String.format("%6.3e", at(i)) + "  ");
			}
		}
		ret = ret.concat("]");
		return ret;
	}
}
