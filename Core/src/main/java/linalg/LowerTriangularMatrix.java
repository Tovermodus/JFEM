package linalg;

import basic.PerformanceArguments;

import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class LowerTriangularMatrix
	implements MutableSolvable
{
	DenseMatrix underlying;
	
	public LowerTriangularMatrix(final int rows, final int cols)
	{
		if (rows != cols)
			throw new IllegalArgumentException("must be square");
		underlying = new DenseMatrix(rows, cols);
	}
	
	public LowerTriangularMatrix(final Matrix matrix)
	{
		if (matrix.getCols() != matrix.getRows())
			throw new IllegalArgumentException("must be square");
		underlying = new DenseMatrix(matrix.getShape());
		if (matrix.isSparse())
		{
			Stream<Map.Entry<IntCoordinates, Double>> entryStream = matrix
				.getCoordinateEntryList()
				.entrySet()
				.stream();
			if (PerformanceArguments.getInstance().parallelizeThreads) entryStream = entryStream.parallel();
			entryStream.filter(e -> e.getKey()
			                         .get(0) >= e.getKey()
			                                     .get(1))
			           .forEach(e -> underlying.add(e.getValue(), e.getKey()));
		} else
		{
			IntStream stream = IntStream.range(0, getRows());
			if (PerformanceArguments.getInstance().parallelizeThreads && matrix.size() > 10000) stream =
				stream.parallel();
			stream.forEach(i ->
			               {
				               for (int j = 0; j <= i; j++)
					               underlying.add(matrix.at(i, j), i, j);
			               });
		}
	}
	
	@Override
	public String toString()
	{
		return printFormatted();
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (!(obj instanceof Matrix)) return false;
		return almostEqual((Tensor) obj);
	}
	
	@Override
	public int hashCode()
	{
		return underlying.hashCode();
	}
	
	@Override
	public double at(final int... coordinates)
	{
		if (coordinates[0] < coordinates[1])
			return 0;
		else return underlying.at(coordinates);
	}
	
	@Override
	public DenseMatrix add(final Tensor other)
	{
		return underlying.add(other);
	}
	
	public LowerTriangularMatrix add(final LowerTriangularMatrix other)
	{
		return new LowerTriangularMatrix(underlying.add(other.underlying));
	}
	
	@Override
	public LowerTriangularMatrix mul(final double scalar)
	{
		return new LowerTriangularMatrix(underlying.mul(scalar));
	}
	
	@Override
	public List<Vector> unfoldDimension(final int dimension)
	{
		return underlying.unfoldDimension(dimension);
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return 0;
	}
	
	@Override
	public boolean isSparse()
	{
		return false;
	}
	
	@Override
	public IntCoordinates getShape()
	{
		return underlying.getShape();
	}
	
	@Override
	public UpperTriangularMatrix transpose()
	{
		return new UpperTriangularMatrix(underlying.transpose());
	}
	
	@Override
	public DenseVector mvMul(final Vector vector)
	{
		return underlying.mvMul(vector);
	}
	
	@Override
	public DenseMatrix mmMul(final Matrix matrix)
	{
		return underlying.mmMul(matrix);
	}
	
	@Override
	public void deleteRow(final int row)
	{
		underlying.deleteRow(row);
	}
	
	@Override
	public void deleteColumn(final int column)
	{
		underlying.deleteColumn(column);
	}
	
	@Override
	public void set(final double value, final int... coordinates)
	{
		if (coordinates[0] >= coordinates[1])
			underlying.set(value, coordinates);
		else throw new IllegalArgumentException("not in lower triangle");
	}
	
	@Override
	public void add(final double value, final int... coordinates)
	{
		if (coordinates[0] >= coordinates[1])
			underlying.add(value, coordinates);
		else throw new IllegalArgumentException("not in lower triangle");
	}
	
	@Override
	public void addInPlace(final Tensor other)
	{
		if (!(other instanceof LowerTriangularMatrix))
			throw new IllegalArgumentException("other needs to be lower triangular");
		underlying.addInPlace(((LowerTriangularMatrix) other).underlying);
	}
	
	@Override
	public void mulInPlace(final double scalar)
	{
		underlying.mulInPlace(scalar);
	}
	
	@Override
	public DenseVector solve(final Vector rhs)
	{
		final DenseVector sol = new DenseVector(rhs.getLength());
		for (int i = 0; i < rhs.size(); i++)
		{
			final double diag = underlying.at(i, i);
			if (diag == 0)
				throw new IllegalStateException("Singular Matrix");
			sol.set(rhs.at(i) / diag, i);
			for (int j = 0; j < i; j++)
			{
				sol.add(-sol.at(j) / diag * underlying.at(i, j), i);
			}
		}
		return sol;
	}
	
	@Override
	public LowerTriangularMatrix inverse()
	{
		return new LowerTriangularMatrix(underlying.inverse());
	}
	
	@Override
	public int getRows()
	{
		return underlying.getRows();
	}
	
	@Override
	public int getCols()
	{
		return underlying.getCols();
	}
	
	@Override
	public Vector solveSymm(final Vector rhs)
	{
		throw new IllegalStateException("is not symm");
	}
}
