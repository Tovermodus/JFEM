package linalg;

import basic.DoubleCompare;
import basic.PerformanceArguments;
import com.google.common.collect.ImmutableMap;

import java.util.*;

public class SparseMatrix implements MutableMatrix, DirectlySolvable, Decomposable
{
	private final int rows;
	private final int cols;
	volatile double[] sparseValues;
	volatile int[] sparseXs;
	volatile int[] sparseYs;
	volatile int sparseEntries;
	
	public SparseMatrix(int rows, int cols)
	{
		this.rows = rows;
		this.cols = cols;
		sparseValues = new double[1];
		sparseYs = new int[1];
		sparseXs = new int[1];
		this.sparseEntries = 1;
	}
	
	public SparseMatrix(SparseMatrix m)
	{
		this.rows = m.rows;
		this.cols = m.cols;
		this.sparseEntries = m.sparseEntries;
		sparseXs = new int[sparseEntries];
		sparseYs = new int[sparseEntries];
		sparseValues = new double[sparseEntries];
		for (int i = 0; i < sparseEntries; i++)
		{
			sparseXs[i] = m.sparseXs[i];
			sparseYs[i] = m.sparseYs[i];
			sparseValues[i] = m.sparseValues[i];
		}
	}
	
	public static SparseMatrix identity(int i)
	{
		SparseMatrix ret = new SparseMatrix(i, i);
		for (int j = 0; j < i; j++)
		{
			ret.add(1, j, j);
		}
		return ret;
	}
	
	
	@Override
	public Vector solve(Vector rhs)
	{
		return new DenseMatrix(this).solve(rhs);
	}
	
	@Override
	public Vector solveSymm(Vector rhs)
	{
		return new DenseMatrix(this).solveSymm(rhs);
	}
	
	@Override
	public double at(int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
		double ret = 0;
		for (int i = 0; i < sparseEntries; i++)
		{
			if (sparseYs[i] == coordinates[0])
				if (sparseXs[i] == coordinates[1])
					ret += sparseValues[i];
		}
		return ret;
	}
	
	public void addSmallMatrixAt(Matrix small, int...coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
			if (coordinates[0]+small.getRows() > getRows())
				throw new IllegalArgumentException("small Matrix too large in y for position");
			if (coordinates[1]+small.getCols() > getCols())
				throw new IllegalArgumentException("small Matrix too large in x for position");
		}
		for(IntCoordinates c:small.getShape().range())
		{
			add(small.at(c), c.get(0) + coordinates[0],
				c.get(1) + coordinates[1]);
		}
		
	}
	public void addSmallMatrixAt(SparseMatrix small, int...coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
			if (coordinates[0]+small.getRows() > getRows())
				throw new IllegalArgumentException("small Matrix too large in y for position");
			if (coordinates[1]+small.getCols() > getCols())
				throw new IllegalArgumentException("small Matrix too large in x for position");
		}
		for(Map.Entry<IntCoordinates, Double> smallEntry: small.getCoordinateEntryList().entrySet())
		{
			add(smallEntry.getValue(), smallEntry.getKey().get(0) + coordinates[0],
				smallEntry.getKey().get(1) + coordinates[1]);
		}
		
	}
	
	private synchronized void resizeSparse()
	{
		TreeMap<IntCoordinates, Double> sparseEntriesMap = new TreeMap<>(getCoordinateEntryList());
		double[] sparseVals = new double[sparseValues.length * 2];
		int[] sparseX = new int[sparseValues.length * 2];
		int[] sparseY = new int[sparseValues.length * 2];
		int i = 0;
		for (Map.Entry<IntCoordinates, Double> entry : sparseEntriesMap.entrySet())
		{
			sparseVals[i] = entry.getValue();
			sparseX[i] = entry.getKey().get(1);
			sparseY[i] = entry.getKey().get(0);
			i++;
		}
		sparseEntries = sparseEntriesMap.size();
		sparseValues = sparseVals;
		sparseYs = sparseY;
		sparseXs = sparseX;
	}
	
	public void deleteRow(int lineCoordinate)
	{
		for (int i = 0; i < sparseEntries; i++)
			if (sparseYs[i] == lineCoordinate)
				sparseValues[i] = 0;
	}
	public void deleteColumn(int lineCoordinate)
	{
		for (int i = 0; i < sparseEntries; i++)
			if (sparseXs[i] == lineCoordinate)
				sparseValues[i] = 0;
	}
	
	@Override
	public synchronized void set(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
		if (sparseEntries >= sparseValues.length - 2)
			resizeSparse();
		sparseValues[sparseEntries] = value - at(coordinates[0], coordinates[1]);
		sparseYs[sparseEntries] = coordinates[0];
		sparseXs[sparseEntries] = coordinates[1];
		sparseEntries++;
	}
	
	@Override
	public synchronized void add(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
			if(coordinates[0] < 0 || coordinates[0] >= getRows())
				throw new IllegalArgumentException("y coordinate out of bounds");
			if(coordinates[1] < 0 || coordinates[1] >= getCols())
				throw new IllegalArgumentException("x coordinate out of bounds");
		}
		if (sparseEntries >= sparseValues.length - 2)
			resizeSparse();
		sparseValues[sparseEntries] = value;
		sparseYs[sparseEntries] = coordinates[0];
		sparseXs[sparseEntries] = coordinates[1];
		sparseEntries++;
	}
	
	@Override
	public void addInPlace(Tensor other)
	{
		
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		if (!(other instanceof SparseMatrix))
		{
			for (IntCoordinates c : other.getShape().range())
				add(other.at(c), c);
		} else
		{
			for (int i = 0; i < ((SparseMatrix) other).sparseEntries; i++)
			{
				add(((SparseMatrix) other).sparseValues[i], ((SparseMatrix) other).sparseYs[i],
					((SparseMatrix) other).sparseXs[i]);
			}
		}
	}
	
	@Override
	public void mulInPlace(double scalar)
	{
		for (int i = 0; i < sparseEntries; i++)
		{
			sparseValues[i] *= scalar;
		}
	}
	
	@Override
	public Matrix add(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		if (!(other instanceof SparseMatrix))
			return ((Matrix) other).add(this);
		return add((SparseMatrix) other);
	}
	
	public SparseMatrix add(SparseMatrix other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		SparseMatrix ret = new SparseMatrix(this);
		for (int i = 0; i < other.sparseEntries; i++)
		{
			ret.add(other.sparseValues[i], other.sparseYs[i],
				other.sparseXs[i]);
		}
		return ret;
	}
	
	public SparseMatrix sub(SparseMatrix other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		SparseMatrix ret = new SparseMatrix(this);
		for (int i = 0; i < other.sparseEntries; i++)
		{
			ret.add(-other.sparseValues[i], other.sparseYs[i],
				other.sparseXs[i]);
		}
		return ret;
	}
	
	@Override
	public ImmutableMap<IntCoordinates, Double> getCoordinateEntryList()
	{
		Map<IntCoordinates, Double> list = new HashMap<>();
		for (int i = 0; i < sparseEntries; i++)
		{
			IntCoordinates coordinates = new IntCoordinates(sparseYs[i], sparseXs[i]);
			double value = sparseValues[i];
			if (list.containsKey(coordinates))
			{
				value += list.get(coordinates);
			}
			list.put(coordinates, value);
		}
		return ImmutableMap.copyOf(list);
	}
	
	@Override
	public SparseMatrix mul(double scalar)
	{
		SparseMatrix ret = new SparseMatrix(this);
		for (int i = 0; i < sparseEntries; i++)
		{
			ret.sparseValues[i] *= scalar;
		}
		return ret;
	}
	
	@Override
	public List<Vector> unfoldDimension(int dimension)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (dimension > 1)
				throw new IllegalArgumentException("Matrix is two dimensional");
		List<Vector> ret = new ArrayList<>(getShape().get(dimension));
		for (int i = 0; i < getShape().get(dimension); i++)
			ret.add(new DenseVector(getShape().get(1 - dimension)));
		if (dimension == 0)
			for (int i = 0; i < sparseEntries; i++)
				((MutableVector) ret.get(sparseYs[i])).add(sparseValues[i], sparseXs[i]);
		if (dimension == 1)
			for (int i = 0; i < sparseEntries; i++)
				((MutableVector) ret.get(sparseXs[i])).add(sparseValues[i], sparseYs[i]);
		
		return ret;
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return sparseEntries;
	}
	
	@Override
	public boolean isSparse()
	{
		return true;
	}
	
	@Override
	public IntCoordinates getShape()
	{
		return new IntCoordinates(rows, cols);
	}
	
	@Override
	public SparseMatrix transpose()
	{
		SparseMatrix ret = new SparseMatrix(cols, rows);
		ret.sparseEntries = sparseEntries;
		ret.sparseValues = sparseValues.clone();
		ret.sparseXs = sparseYs.clone();
		ret.sparseYs = sparseXs.clone();
		return ret;
	}
	
	@Override
	public Vector mvMul(Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (vector.getLength()))
				throw new IllegalArgumentException("Incompatible sizes");
		DenseVector ret = new DenseVector(getRows());
		for (int i = 0; i < sparseEntries; i++)
			ret.add(sparseValues[i] * vector.at(sparseXs[i]), sparseYs[i]);
		return ret;
		
	}
	
	@Override
	public Vector tvMul(Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (vector.getLength()))
				throw new IllegalArgumentException("Incompatible sizes");
		DenseVector ret = new DenseVector(getCols());
		for (int i = 0; i < sparseEntries; i++)
			ret.add(sparseValues[i] * vector.at(sparseYs[i]), sparseXs[i]);
		return ret;
		
	}
	
	@Override
	public DirectlySolvable mmMul(Matrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (matrix.getRows()))
				throw new IllegalArgumentException("Incompatible sizes");
		if(! (matrix instanceof  SparseMatrix))
			return (new DenseMatrix(this)).mmMul(matrix);
		SparseMatrix ret = new SparseMatrix(getRows(), matrix.getCols());
		for (int i = 0; i < sparseEntries; i++)
			for (int j = 0; j < ((SparseMatrix) matrix).sparseEntries; j++)
			{
				if (sparseXs[i] == ((SparseMatrix) matrix).sparseYs[j])
					ret.add(sparseValues[i] * ((SparseMatrix) matrix).sparseValues[j], sparseYs[i],
						((SparseMatrix) matrix).sparseXs[j]);
			}
		return ret;
	}
	
	@Override
	public SparseMatrix getLowerTriangleMatrix()
	{
		SparseMatrix mat = new SparseMatrix(getRows(), getCols());
		for (int i = 0; i < sparseEntries; i++)
		{
			if (sparseYs[i] > sparseXs[i])
				mat.add(sparseValues[i], sparseYs[i], sparseXs[i]);
		}
		return mat;
	}
	
	@Override
	public SparseMatrix getUpperTriangleMatrix()
	{
		SparseMatrix mat = new SparseMatrix(getRows(), getCols());
		for (int i = 0; i < sparseEntries; i++)
		{
			if (sparseYs[i] < sparseXs[i])
				mat.add(sparseValues[i], sparseYs[i], sparseXs[i]);
		}
		return mat;
	}
	
	@Override
	public boolean almostEqual(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Tensors are of different size");
		}
		if (!(other instanceof SparseMatrix))
			return other.almostEqual(this);
		ImmutableMap<IntCoordinates, Double> myValues = getCoordinateEntryList();
		ImmutableMap<IntCoordinates, Double> otherValues = other.getCoordinateEntryList();
		double absmax = absMaxElement() + other.absMaxElement();
		for (Map.Entry<IntCoordinates, Double> entry : myValues.entrySet())
		{
			if (otherValues.containsKey(entry.getKey()))
			{
				if (!DoubleCompare.almostEqualAfterOps(entry.getValue(), otherValues.get(entry.getKey()),
					absmax*10, sparseEntries))
				{
					System.out.println(entry.getValue() + " != " + otherValues.get(entry.getKey()) +" difference: " + Math.abs(entry.getValue() - otherValues.get(entry.getKey())) + " at " + entry.getKey() + " out of " + getShape());
					return false;
				}
			} else
			{
				System.out.println("Value at" + entry.getKey() + " was not found in second matrix");
				return false;
			}
		}
		if(otherValues.size() != myValues.size())
			System.out.println("matrices have different numbers of values");
		return otherValues.size() == myValues.size();
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if (!(obj instanceof Matrix))
			return false;
		return almostEqual((Tensor) obj);
	}
	
	@Override
	public String toString()
	{
		return printFormatted();
	}
	
	@Override
	public SparseMatrix getDiagonalMatrix()
	{
		SparseMatrix mat = new SparseMatrix(getRows(), getCols());
		for (int i = 0; i < sparseEntries; i++)
		{
			if (sparseYs[i] == sparseXs[i])
				mat.add(sparseValues[i], sparseYs[i], sparseXs[i]);
		}
		return mat;
	}
}
