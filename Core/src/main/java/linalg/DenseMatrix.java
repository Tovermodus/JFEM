package linalg;

import basic.DoubleCompare;
import basic.PerformanceArguments;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class DenseMatrix
	implements MutableSolvable, Decomposable, Serializable
{
	protected volatile double[][] entries;
	private transient org.ujmp.core.Matrix ujmpmat = null;
	
	public DenseMatrix(final int i, final int j)
	{
		entries = new double[i][j];
	}
	
	public DenseMatrix(final Matrix matrix)
	{
		entries = new double[matrix.getRows()][matrix.getCols()];
		if (matrix.isSparse())
		{
			Stream<Map.Entry<IntCoordinates, Double>> entryStream = matrix
				.getCoordinateEntryList()
				.entrySet()
				.stream();
			if (PerformanceArguments.getInstance().parallelizeThreads) entryStream = entryStream.parallel();
			entryStream.forEach(e -> entries[e.getKey()
			                                  .get(0)][e.getKey()
			                                            .get(1)] = e.getValue());
		} else
		{
			IntStream stream = IntStream.range(0, getRows());
			if (PerformanceArguments.getInstance().parallelizeThreads && matrix.size() > 10000) stream =
				stream.parallel();
			stream.forEach(i ->
			               {
				               for (int j = 0; j < getCols(); j++)
					               entries[i][j] = matrix.at(i, j);
			               });
		}
	}
	
	public DenseMatrix(final double[][] matrix)
	{
		entries = matrix;
	}
	
	public DenseMatrix(final DenseMatrix d, final boolean wrap)
	{
		if (wrap) this.entries = d.entries;
		else
		{
			entries = new double[d.getRows()][d.getCols()];
			for (int i = 0; i < getRows(); i++)
				for (int j = 0; j < getCols(); j++)
					entries[i][j] = d.at(i, j);
		}
	}
	
	public DenseMatrix(final IntCoordinates rowCols)
	{
		this(rowCols.get(0), rowCols.get(1));
	}
	
	public static DenseMatrix squareMatrixFromValues(final double... values)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (Math.pow((int) Math.sqrt(values.length), 2) != values.length)
				throw new UnsupportedOperationException("nope");
		final DenseMatrix ret = new DenseMatrix((int) Math.sqrt(values.length), (int) Math.sqrt(values.length));
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(values[i * ret.getRows() + j], i, j);
		return ret;
	}
	
	public static DenseMatrix identity(final int n)
	{
		final DenseMatrix ret = new DenseMatrix(n, n);
		for (int i = 0; i < n; i++)
			ret.add(1, i, i);
		return ret;
	}
	
	private static DenseMatrix fromUJMPMatrix(final org.ujmp.core.Matrix matrix)
	{
		final DenseMatrix ret = new DenseMatrix((int) matrix.getRowCount(), (int) matrix.getColumnCount());
		
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(matrix.getAsDouble(i, j), i, j);
		return ret;
	}
	
	private static SparseMatrix fromUJMPMatrixSparse(final org.ujmp.core.Matrix matrix)
	{
		final SparseMatrix ret = new SparseMatrix((int) matrix.getRowCount(), (int) matrix.getColumnCount());
		
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				if (Math.abs(
					matrix.getAsDouble(i, j)) > PerformanceArguments.getInstance().doubleTolerance)
					ret.add(matrix.getAsDouble(i, j), i, j);
		return ret;
	}
	
	@Override
	public double frobeniusInner(final Matrix other)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		double ret = 0;
		for (int i = 0; i < getRows(); i++)
			for (int j = 0; j < getCols(); j++)
				ret += at(i, j) * other.at(i, j);
		return ret;
	}
	
	@Override
	public double at(final int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2) throw new IllegalArgumentException("Wrong number of coordinates");
		return entries[coordinates[0]][coordinates[1]];
	}
	
	@Override
	public void set(final double value, final int... coordinates)
	{
		//if (PerformanceArguments.getInstance().executeChecks)
		if (coordinates.length != 2) throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]][coordinates[1]] = value;
	}
	
	@Override
	public void add(final double value, final int... coordinates)
	{
		//if (PerformanceArguments.getInstance().executeChecks)
		if (coordinates.length != 2) throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]][coordinates[1]] += value;
	}
	
	@Override
	public int getRows()
	{
		return entries.length;
	}
	
	@Override
	public int getCols()
	{
		return entries[0].length;
	}
	
	@Override
	public void addInPlace(final Tensor other)
	{
		if (other instanceof Matrix && other.isSparse())
			for (final var entry : other.getCoordinateEntryList()
			                            .entrySet())
				add(entry.getValue(), entry.getKey());
		else
			for (int i = 0; i < getRows(); i++)
				for (int j = 0; j < getCols(); j++)
					entries[i][j] += other.at(i, j);
	}
	
	@Override
	public void subInPlace(final Tensor other)
	{
		if (other instanceof Matrix && other.isSparse())
			for (final var entry : other.getCoordinateEntryList()
			                            .entrySet())
				add(-entry.getValue(), entry.getKey());
		else
			for (int i = 0; i < getRows(); i++)
				for (int j = 0; j < getCols(); j++)
					entries[i][j] -= other.at(i, j);
	}
	
	@Override
	public void mulInPlace(final double scalar)
	{
		for (int i = 0; i < getRows(); i++)
			for (int j = 0; j < getCols(); j++)
				entries[i][j] *= scalar;
	}
	
	@Override
	public DenseMatrix sub(final Tensor other)
	{
		
		if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		if (other instanceof Matrix && other.isSparse())
		{
			if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
			final DenseMatrix ret = new DenseMatrix(this);
			other.getCoordinateEntryList()
			     .forEach((k, v) -> ret.add(-v, k));
			return ret;
		} else
		{
			final DenseMatrix ret = new DenseMatrix(this);
			for (int i = 0; i < ret.getRows(); i++)
				for (int j = 0; j < ret.getCols(); j++)
					ret.add(-other.at(i, j), i, j);
			return ret;
		}
	}
	
	@Override
	public DenseMatrix slice(final IntCoordinates start, final IntCoordinates end)
	{
		final DenseMatrix ret = new DenseMatrix(end.sub(start));
		for (int i = start.get(0); i < end.get(0); i++)
			for (int j = start.get(1); j < end.get(1); j++)
				ret.entries[i - start.get(0)][j - start.get(1)] = entries[i][j];
		return ret;
	}
	
	@Override
	public DenseMatrix add(final Tensor other)
	{
		
		if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		if (other instanceof Matrix && other.isSparse())
		{
			if (PerformanceArguments.getInstance().executeChecks) if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
			final DenseMatrix ret = new DenseMatrix(this);
			other.getCoordinateEntryList()
			     .forEach((k, v) -> ret.add(v, k));
			return ret;
		} else
		{
			final DenseMatrix ret = new DenseMatrix(this);
			for (int i = 0; i < ret.getRows(); i++)
				for (int j = 0; j < ret.getCols(); j++)
					ret.add(other.at(i, j), i, j);
			return ret;
		}
	}
	
	@Override
	public DenseMatrix mul(final double scalar)
	{
		
		final DenseMatrix ret = new DenseMatrix(entries.length, entries[0].length);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(scalar * at(i, j), i, j);
		return ret;
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
				final DenseVector vec = new DenseVector(getCols());
				for (int j = 0; j < getCols(); j++)
					vec.add(at(i, j), j);
				ret.add(vec);
			}
		} else
		{
			for (int i = 0; i < getCols(); i++)
			{
				final DenseVector vec = new DenseVector(getRows());
				for (int j = 0; j < getRows(); j++)
					vec.add(at(j, i), j);
				ret.add(vec);
			}
		}
		return ret;
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
		return new IntCoordinates(entries.length, entries[0].length);
	}
	
	@Override
	public long size()
	{
		return (long) entries.length * entries[0].length;
	}
	
	public boolean almostEqual(final SparseMatrix other)
	{
		return almostEqual(new DenseMatrix(other));
	}
	
	@Override
	public DenseMatrix transpose()
	{
		final DenseMatrix ret = new DenseMatrix(entries[0].length, entries.length);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(at(j, i), i, j);
		return ret;
	}
	
	@Override
	public DenseVector mvMul(final Vector v)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (v.getLength())) throw new IllegalArgumentException("Incompatible sizes");
		final int cols = getCols();
		final DenseVector ret = new DenseVector(getRows());
		final DenseVector vector = new DenseVector(v);
		for (int i = 0; i < getRows(); i++)
		{
			double ret_i = 0;
			final double[] row = entries[i];
			for (int k = 0; k < cols; k++)
				ret_i += row[k] * vector.entries[k];
			ret.add(ret_i, i);
		}
		return ret;
	}
	
	@Override
	public DenseVector tvMul(final Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != (vector.getLength())) throw new IllegalArgumentException("Incompatible sizes");
		final DenseVector ret = new DenseVector(getCols());
		for (int k = 0; k < getRows(); k++)
			for (int i = 0; i < getCols(); i++)
				ret.add(at(k, i) * vector.at(k), i);
		return ret;
	}
	
	public DenseMatrix mmMul(final DenseMatrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (matrix.getRows())) throw new IllegalArgumentException("Incompatible sizes");
		final DenseMatrix ret = new DenseMatrix(getRows(), matrix.getCols());
		IntStream str = IntStream.range(0, getRows());
		if (PerformanceArguments.getInstance().parallelizeThreads)
			str = str.parallel();
		final int cols = getCols();
		final int otherCols = matrix.getCols();
		str.forEach(i ->
		            {
			            final double[] retRow = ret.entries[i];
			            for (int k = 0; k < cols; k++)
			            {
				            final double aik = at(i, k);
				            final double[] matRow = matrix.entries[k];
				            for (int j = 0; j < otherCols; j++)
					            retRow[j] += aik * matRow[j];
			            }
		            });
		return ret;
	}
	
	private void overrideBy(final DenseMatrix other)
	{
		this.entries = other.entries;
		this.ujmpmat = other.ujmpmat;
	}
	
	@Override
	public DenseMatrix mmMul(final Matrix matrix)
	{
		if (!matrix.isSparse())
			return mmMul(new DenseMatrix(matrix));
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (matrix.getRows())) throw new IllegalArgumentException("Incompatible sizes");
		final IntList[] columnCoordinates = new IntList[matrix.getCols()];
		final DoubleList[] columnEntries = new DoubleList[matrix.getCols()];
		groupCoordinateEntriesByCol(matrix, columnCoordinates, columnEntries);
		final DenseMatrix ret = new DenseMatrix(getRows(), matrix.getCols());
		final int rows = getRows();
		for (int i = 0; i < rows; i++)
		{
			final double[] row = entries[i];
			final int finalI = i;
			IntStream str = IntStream.range(0, matrix.getCols());
			if (PerformanceArguments.getInstance().parallelizeThreads)
				str = str.parallel();
			str.forEach(columnIndex ->
			            {
				            if (columnCoordinates[columnIndex] != null)
				            {
					            final int[] columnCoords
						            = columnCoordinates[columnIndex].toIntArray();
					            final double[] columnEntriees
						            = columnEntries[columnIndex].toDoubleArray();
					            double contraction = 0;
					            for (int j = 0; j < columnCoords.length; j++)
						            contraction
							            += row[columnCoords[j]] * columnEntriees[j];
					            ret.add(contraction, finalI, columnIndex);
				            }
			            });
		}
		return ret;
	}
	
	public static void groupCoordinateEntriesByCol(final Matrix matrix,
	                                               final IntList[] columnCoordinates,
	                                               final DoubleList[] columnEntries)
	{
		matrix.getCoordinateEntryList()
		      .forEach((intCoordinates, value) ->
		               {
			               final int col = intCoordinates.get(1);
			               if (columnCoordinates[col] == null)
			               {
				               columnCoordinates[col] = new IntArrayList();
				               columnEntries[col] = new DoubleArrayList();
			               }
			               columnCoordinates[col].add(intCoordinates.get(0));
			               columnEntries[col].add(value);
		               });
	}
	
	@Override
	public DenseMatrix getStrictlyLowerTriangleMatrix()
	{
		if (PerformanceArguments.getInstance().executeChecks) if (getRows() != getCols())
			throw new UnsupportedOperationException("cant get triangle from rectangle");
		final DenseMatrix ret = new DenseMatrix(getRows(), getCols());
		for (int i = 0; i < getRows(); i++)
			for (int j = 0; j < i; j++)
				ret.add(at(i, j), i, j);
		return ret;
	}
	
	@Override
	public DenseMatrix getStrictlyUpperTriangleMatrix()
	{
		if (PerformanceArguments.getInstance().executeChecks) if (getRows() != getCols())
			throw new UnsupportedOperationException("cant get triangle from rectangle");
		final DenseMatrix ret = new DenseMatrix(getRows(), getCols());
		for (int i = 0; i < getRows(); i++)
			for (int j = i + 1; j < getCols(); j++)
				ret.add(at(i, j), i, j);
		return ret;
	}
	
	@Override
	public DenseMatrix getDiagonalMatrix()
	{
		if (PerformanceArguments.getInstance().executeChecks) if (getRows() != getCols())
			throw new UnsupportedOperationException("cant get triangle from rectangle");
		final DenseMatrix ret = new DenseMatrix(getRows(), getCols());
		for (int i = 0; i < getRows(); i++)
			ret.add(at(i, i), i, i);
		return ret;
	}
	
	private org.ujmp.core.Matrix toUJMPMatrix()
	{
		if (ujmpmat == null)
		{
			ujmpmat = org.ujmp.core.DenseMatrix.Factory.zeros(getRows(), getCols());
			for (int i = 0; i < getRows(); i++)
				for (int j = 0; j < getCols(); j++)
					ujmpmat.setAsDouble(at(i, j), i, j);
		}
		return ujmpmat;
	}
	
	public SparseMatrix getCholeskyL()
	{
		final org.ujmp.core.Matrix c = org.ujmp.core.Matrix.chol.calc(toUJMPMatrix());
		return fromUJMPMatrixSparse(c);
	}
	
	@Override
	public DenseMatrix inverse()
	{
		return fromUJMPMatrix(toUJMPMatrix().inv());
	}
	
	@Override
	public DenseVector solve(final Vector rhs)
	{
		return DenseVector.fromUJMPVector(toUJMPMatrix().solve(new DenseVector(rhs).toUJMPVector()));
	}
	
	@Override
	public Vector solveSymm(final Vector rhs)
	{
		return DenseVector.fromUJMPVector(toUJMPMatrix().solveSymm(new DenseVector(rhs).toUJMPVector()));
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
		return Arrays
			.stream(entries)
			.mapToInt(e -> Arrays.stream(e)
			                     .mapToInt(DoubleCompare::doubleHash)
			                     .sum())
			.sum();
	}
	
	@Override
	public String toString()
	{
		return printFormatted();
	}
	
	@Override
	public void deleteRow(final int row)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (0 > row || row > getRows()) throw new UnsupportedOperationException("row out of bounds");
		for (int i = 0; i < getCols(); i++)
			entries[row][i] = 0;
	}
	
	@Override
	public void deleteColumn(final int column)
	{
		
		if (PerformanceArguments.getInstance().executeChecks) if (0 > column || column > getCols())
			throw new UnsupportedOperationException("column out of bounds");
		for (int i = 0; i < getRows(); i++)
			entries[i][column] = 0;
	}
	
	public RealMatrix toApacheMatrix()
	{
		return new Array2DRowRealMatrix(entries);
	}
}
