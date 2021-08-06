package linalg;

import basic.DoubleCompare;
import basic.PerformanceArguments;
import basic.PlotWindow;
import com.google.common.primitives.Ints;
import org.ujmp.core.doublematrix.calculation.general.decomposition.Chol;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

public class DenseMatrix implements MutableMatrix, Decomposable, DirectlySolvable
{
	protected volatile double[][] entries;
	private org.ujmp.core.Matrix ujmpmat = null;
	
	public DenseMatrix(int i, int j)
	{
		entries = new double[i][j];
	}
	
	public DenseMatrix(Matrix matrix)
	{
		entries = new double[matrix.getRows()][matrix.getCols()];
		IntStream stream = IntStream.range(0,getRows());
		if(PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(i ->
		{
			for (int j = 0; j < getCols(); j++)
				entries[i][j] = matrix.at(i, j);
		});
	}
	
	public DenseMatrix(double [][] matrix)
	{
		entries = matrix;
	}
	
	public DenseMatrix(DenseMatrix d, boolean wrap)
	{
		if(wrap)
			this.entries = d.entries;
		else
		{
			entries = new double[d.getRows()][d.getCols()];
			for (int i = 0; i < getRows(); i++)
				for (int j = 0; j < getCols(); j++)
					entries[i][j] = d.at(i, j);
		}
	}
	
	@Override
	public double frobeniusInner(Matrix other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		double ret = 0;
		for (int i = 0; i < getRows(); i++)
			for (int j = 0; j < getCols(); j++)
				ret += at(i, j) * other.at(i, j);
		return ret;
	}
	
	public static DenseMatrix squareMatrixFromValues(double... values)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (Math.pow((int) Math.sqrt(values.length), 2) != values.length)
				throw new UnsupportedOperationException("nope");
		DenseMatrix ret = new DenseMatrix((int) Math.sqrt(values.length), (int) Math.sqrt(values.length));
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(values[i * ret.getRows() + j], i, j);
		return ret;
	}
	
	public static DenseMatrix identity(int n)
	{
		DenseMatrix ret = new DenseMatrix(n, n);
		for (int i = 0; i < n; i++)
			ret.add(1, i, i);
		return ret;
	}
	
	@Override
	public double at(int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
		return entries[coordinates[0]][coordinates[1]];
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]][coordinates[1]] = value;
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]][coordinates[1]] += value;
	}
	
	@Override
	public void addInPlace(Tensor other)
	{
		for (int i = 0; i < getRows(); i++)
			for (int j = 0; j < getCols(); j++)
				entries[i][j] += other.at(i, j);
	}
	
	@Override
	public void subInPlace(Tensor other)
	{
		for (int i = 0; i < getRows(); i++)
			for (int j = 0; j < getCols(); j++)
				entries[i][j] -= other.at(i, j);
	}
	
	@Override
	public void mulInPlace(double scalar)
	{
		for (int i = 0; i < getRows(); i++)
			for (int j = 0; j < getCols(); j++)
				entries[i][j] *= scalar;
	}
	
	@Override
	public DenseMatrix add(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		DenseMatrix ret = new DenseMatrix(this);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.add(other.at(i, j), i, j);
		return ret;
	}
	
	@Override
	public DenseMatrix sub(Tensor other)
	{
		
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		DenseMatrix ret = new DenseMatrix(this);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.add(-other.at(i, j), i, j);
		return ret;
	}
	
	@Override
	public DenseMatrix mul(double scalar)
	{
		
		DenseMatrix ret = new DenseMatrix(entries.length, entries[0].length);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(scalar * at(i, j), i, j);
		return ret;
	}
	
	@Override
	public List<Vector> unfoldDimension(int dimension)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (dimension > 1)
				throw new IllegalArgumentException("Matrix is two dimensional");
		List<Vector> ret = new ArrayList<>();
		if (dimension == 0)
		{
			for (int i = 0; i < getRows(); i++)
			{
				DenseVector vec = new DenseVector(getCols());
				for (int j = 0; j < getCols(); j++)
					vec.add(at(i, j), j);
				ret.add(vec);
			}
		} else
		{
			for (int i = 0; i < getCols(); i++)
			{
				DenseVector vec = new DenseVector(getRows());
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
		return entries.length * entries[0].length;
	}
	
	
	@Override
	public DenseMatrix transpose()
	{
		DenseMatrix ret = new DenseMatrix(entries[0].length, entries.length);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(at(j, i), i, j);
		return ret;
	}
	
	@Override
	public DenseVector mvMul(Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (vector.getLength()))
				throw new IllegalArgumentException("Incompatible sizes");
		DenseVector ret = new DenseVector(getRows());
		for (int i = 0; i < getRows(); i++)
			for (int k = 0; k < getCols(); k++)
				ret.add(at(i, k) * vector.at(k), i);
		return ret;
	}
	
	@Override
	public DenseVector tvMul(Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != (vector.getLength()))
				throw new IllegalArgumentException("Incompatible sizes");
		DenseVector ret = new DenseVector(getCols());
		for (int i = 0; i < getCols(); i++)
			for (int k = 0; k < getRows(); k++)
				ret.add(at(k, i) * vector.at(k), i);
		return ret;
	}
	
	@Override
	public DenseMatrix mmMul(Matrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (matrix.getRows()))
				throw new IllegalArgumentException("Incompatible sizes");
		DenseMatrix ret = new DenseMatrix(getRows(), matrix.getCols());
		for (int i = 0; i < getRows(); i++)
			for (int j = 0; j < matrix.getCols(); j++)
				for (int k = 0; k < getCols(); k++)
					ret.add(at(i, k) * matrix.at(k, j), i, j);
		return ret;
	}
	
	
	@Override
	public DenseMatrix getLowerTriangleMatrix()
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != getCols())
				throw new UnsupportedOperationException("cant get triangle from rectangle");
		DenseMatrix ret = new DenseMatrix(getRows(), getCols());
		for (int i = 0; i < getRows(); i++)
			for (int j = 0; j < i; j++)
				ret.add(at(i, j), i, j);
		return ret;
	}
	
	@Override
	public DenseMatrix getUpperTriangleMatrix()
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != getCols())
				throw new UnsupportedOperationException("cant get triangle from rectangle");
		DenseMatrix ret = new DenseMatrix(getRows(), getCols());
		for (int i = 0; i < getRows(); i++)
			for (int j = i + 1; j < getCols(); j++)
				ret.add(at(i, j), i, j);
		return ret;
	}
	
	@Override
	public DenseMatrix getDiagonalMatrix()
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != getCols())
				throw new UnsupportedOperationException("cant get triangle from rectangle");
		DenseMatrix ret = new DenseMatrix(getRows(), getCols());
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
		org.ujmp.core.Matrix c = org.ujmp.core.Matrix.chol.calc(toUJMPMatrix());
		return fromUJMPMatrixSparse(c);
	}
	
	
	private static DenseMatrix fromUJMPMatrix(org.ujmp.core.Matrix matrix)
	{
		DenseMatrix ret = new DenseMatrix((int) matrix.getRowCount(), (int) matrix.getColumnCount());
		
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(matrix.getAsDouble(i, j), i, j);
		return ret;
	}
	private static SparseMatrix fromUJMPMatrixSparse(org.ujmp.core.Matrix matrix)
	{
		SparseMatrix ret = new SparseMatrix((int) matrix.getRowCount(), (int) matrix.getColumnCount());
		
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				if(Math.abs(matrix.getAsDouble(i,j)) > PerformanceArguments.getInstance().doubleTolerance)
					ret.add(matrix.getAsDouble(i, j), i, j);
		return ret;
	}
	
	@Override
	public DenseMatrix inverse()
	{
		return fromUJMPMatrix(toUJMPMatrix().inv());
	}
	
	@Override
	public DenseVector solve(Vector rhs)
	{
		return DenseVector.fromUJMPVector(toUJMPMatrix().solve(new DenseVector(rhs).toUJMPVector()));
	}
	
	@Override
	public Vector solveSymm(Vector rhs)
	{
		return DenseVector.fromUJMPVector(toUJMPMatrix().solveSymm(new DenseVector(rhs).toUJMPVector()));
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if (!(obj instanceof Matrix))
			return false;
		return almostEqual((Tensor) obj);
	}
	@Override
	public int hashCode()
	{
		return Arrays.stream(entries).mapToInt(e -> Arrays.stream(e).mapToInt(DoubleCompare::doubleHash).sum()).sum();
	}
	
	@Override
	public String toString()
	{
		return printFormatted();
	}
	
	@Override
	public void deleteRow(int row)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (0 <= row && row <= getRows())
				throw new UnsupportedOperationException("row out of bounds");
		for(int i = 0; i < getCols(); i++)
			entries[row][i] = 0;
	}
	
	@Override
	public void deleteColumn(int column)
	{
		
		if (PerformanceArguments.getInstance().executeChecks)
			if (0 <= column && column <= getCols())
				throw new UnsupportedOperationException("column out of bounds");
		for(int i = 0; i < getRows(); i++)
			entries[i][column] = 0;
	}
}
