package linalg;

import java.util.List;

public class LowerTriangularSparseMatrix
	implements MutableSolvable
{
	SparseMatrix underlying;
	SparseMvMul mvm;
	
	public LowerTriangularSparseMatrix(final int rows, final int cols)
	{
		if (rows != cols)
			throw new IllegalArgumentException("must be square");
		underlying = new SparseMatrix(rows, cols);
	}
	
	public LowerTriangularSparseMatrix(final Matrix matrix)
	{
		if (matrix.getCols() != matrix.getRows())
			throw new IllegalArgumentException("must be square");
		underlying = new SparseMatrix(matrix.getRows(), matrix.getCols(), matrix.getSparseEntryCount());
		if (matrix.isSparse())
		{
			final SparseMatrix.ElementOperation op = new SparseMatrix.ElementOperation()
			{
				@Override
				public void operation(final int column, final int row, final double value)
				{
					if (column >= row)
						underlying.add(value, column, row);
				}
			};
			matrix.forEachElement(op);
		} else
			throw new IllegalArgumentException("Use Triangular Dense Matrix");
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
	public SparseMatrix add(final Tensor other)
	{
		return underlying.add(other);
	}
	
	public LowerTriangularSparseMatrix add(final LowerTriangularSparseMatrix other)
	{
		return new LowerTriangularSparseMatrix(underlying.add(other.underlying));
	}
	
	@Override
	public LowerTriangularSparseMatrix mul(final double scalar)
	{
		return new LowerTriangularSparseMatrix(underlying.mul(scalar));
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
		return true;
	}
	
	@Override
	public IntCoordinates getShape()
	{
		return underlying.getShape();
	}
	
	@Override
	public UpperTriangularSparseMatrix transpose()
	{
		return new UpperTriangularSparseMatrix(underlying.transpose());
	}
	
	@Override
	public DenseVector mvMul(final Vector vector)
	{
		if (mvm == null)
			mvm = new SparseMvMul(underlying);
		return mvm.mvMul(vector);
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
		mvm = null;
		if (coordinates[0] >= coordinates[1])
			underlying.set(value, coordinates);
		else throw new IllegalArgumentException("not in lower triangle");
	}
	
	@Override
	public void add(final double value, final int... coordinates)
	{
		mvm = null;
		if (coordinates[0] >= coordinates[1])
			underlying.add(value, coordinates);
		else throw new IllegalArgumentException("not in lower triangle");
	}
	
	@Override
	public void addInPlace(final Tensor other)
	{
		mvm = null;
		if (!(other instanceof UpperTriangularSparseMatrix))
			throw new IllegalArgumentException("other needs to be lower triangular");
		underlying.addInPlace(((LowerTriangularSparseMatrix) other).underlying);
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
		if (mvm == null)
		{
			System.out.println("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMmm");
			mvm = new SparseMvMul(underlying);
		}
		for (int i = 0; i < rhs.getLength(); i++)
		{
			backsubstituteRow(rhs, sol, i);
		}
		return sol;
	}
	
	private void backsubstituteRow(final Vector rhs, final DenseVector sol, final int i)
	{
		final var row = mvm.getRow(i);
		final int[] indices = row._1;
		final double[] vals = row._2;
		double diag = 0;
		int diagIndex = -1;
		for (int j = 0; j < indices.length; j++)
			if (indices[j] == i)
			{
				diag = vals[j];
				diagIndex = j;
				break;
			}
		if (diag == 0)
			throw new IllegalStateException("Singular Matrix");
		sol.set(rhs.at(i) / diag, i);
		
		for (int j = 0; j < indices.length; j++)
		{
			if (j == diagIndex)
				continue;
			sol.add(-sol.at(indices[j]) / diag * vals[j], i);
		}
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
