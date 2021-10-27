package linalg;

import basic.DoubleCompare;
import basic.MapKeySelectCollector;
import basic.PerformanceArguments;
import com.google.common.collect.ImmutableMap;
import io.vavr.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SparseMatrix
	implements MutableMatrix, DirectlySolvable, Decomposable, Serializable
{
	private final int rows;
	private final int cols;
	volatile double[] sparseValues;
	volatile int[] sparseXs;
	volatile int[] sparseYs;
	volatile int sparseEntries;
	
	@Override
	public SparseMatrix slice(final IntCoordinates start, final IntCoordinates end)
	{
		final SparseMatrix ret = new SparseMatrix(end.get(0) - start.get(0), end.get(1) - start.get(1));
		for (int i = 0; i < sparseEntries; i++)
		{
			if (!(sparseXs[i] >= end.get(1)
				      || sparseXs[i] < start.get(1)
				      || sparseYs[i] >= end.get(0)
				      || sparseYs[i] < start.get(0)))
				ret.add(sparseValues[i], sparseYs[i] - start.get(0), sparseXs[i] - start.get(1));
		}
		return ret;
	}
	
	public SparseMatrix(final int rows, final int cols)
	{
		this.rows = rows;
		this.cols = cols;
		sparseValues = new double[1];
		sparseYs = new int[1];
		sparseXs = new int[1];
		this.sparseEntries = 1;
	}
	
	public SparseMatrix(final SparseMatrix m)
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
	
	public SparseMatrix(final Matrix m)
	{
		this.rows = m.getRows();
		this.cols = m.getCols();
		sparseValues = new double[1];
		sparseYs = new int[1];
		sparseXs = new int[1];
		this.sparseEntries = 1;
		resetFromEntries(m.getCoordinateEntryList());
	}
	
	public static SparseMatrix identity(final int i)
	{
		final SparseMatrix ret = new SparseMatrix(i, i);
		for (int j = 0; j < i; j++)
		{
			ret.add(1, j, j);
		}
		return ret;
	}
	
	@Override
	public DenseVector solve(final Vector rhs)
	{
		return new DenseMatrix(this).solve(rhs);
	}
	
	@Override
	public Vector solveSymm(final Vector rhs)
	{
		return new DenseMatrix(this).solveSymm(rhs);
	}
	
	@Override
	public double at(final int... coordinates)
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
	
	public void addSmallMatrixAt(final SparseMatrix small, final int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
			if (coordinates[0] + small.getRows() > getRows())
				throw new IllegalArgumentException("small Matrix too large in y for position");
			if (coordinates[1] + small.getCols() > getCols())
				throw new IllegalArgumentException("small Matrix too large in x for position");
		}
		for (final Map.Entry<IntCoordinates, Double> smallEntry : small.getCoordinateEntryList()
		                                                               .entrySet())
		{
			add(smallEntry.getValue(),
			    smallEntry.getKey()
			              .get(0) + coordinates[0],
			    smallEntry.getKey()
			              .get(1) + coordinates[1]);
		}
	}
	
	private synchronized void resizeSparse()
	{
		final TreeMap<IntCoordinates, Double> sparseEntriesMap = new TreeMap<>(getCoordinateEntryList());
		resetFromEntries(sparseEntriesMap);
	}
	
	private void resetFromEntries(final Map<IntCoordinates, Double> sparseEntriesMap)
	{
		final double[] sparseVals = new double[sparseEntriesMap.size() * 2];
		final int[] sparseX = new int[sparseEntriesMap.size() * 2];
		final int[] sparseY = new int[sparseEntriesMap.size() * 2];
		int i = 0;
		for (final Map.Entry<IntCoordinates, Double> entry : sparseEntriesMap.entrySet())
		{
			sparseVals[i] = entry.getValue();
			sparseX[i] = entry.getKey()
			                  .get(1);
			sparseY[i] = entry.getKey()
			                  .get(0);
			i++;
		}
		sparseEntries = sparseEntriesMap.size();
		sparseValues = sparseVals;
		sparseYs = sparseY;
		sparseXs = sparseX;
	}
	
	@Override
	public void deleteRow(final int lineCoordinate)
	{
		for (int i = 0; i < sparseEntries; i++)
			if (sparseYs[i] == lineCoordinate)
				sparseValues[i] = 0;
	}
	
	@Override
	public void deleteColumn(final int lineCoordinate)
	{
		for (int i = 0; i < sparseEntries; i++)
			if (sparseXs[i] == lineCoordinate)
				sparseValues[i] = 0;
	}
	
	@Override
	public void addColumn(final Vector vector, final int column, final int start)
	{
		
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (start + vector.getLength() > getRows())
				throw new IllegalArgumentException("small Vector too large for position");
		}
		for (int i = 0; i < vector.getLength(); i++)
		{
			if (vector.at(i) != 0)
				add(vector.at(i), i + start, column);
		}
	}
	
	@Override
	public synchronized void set(final double value, final int... coordinates)
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
	public synchronized void add(final double value, final int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
			if (coordinates[0] < 0 || coordinates[0] >= getRows())
				throw new IllegalArgumentException("y coordinate out of bounds" + coordinates[0] + " " +
					                                   "with rows: " + getRows());
			if (coordinates[1] < 0 || coordinates[1] >= getCols())
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
	public void addInPlace(final Tensor other)
	{
		
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		if (!(other instanceof SparseMatrix))
		{
			for (final IntCoordinates c : other.getShape()
			                                   .range())
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
	public void mulInPlace(final double scalar)
	{
		for (int i = 0; i < sparseEntries; i++)
		{
			sparseValues[i] *= scalar;
		}
	}
	
	@Override
	public Matrix add(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		if (!(other instanceof SparseMatrix))
			return ((Matrix) other).add(this);
		return add((SparseMatrix) other);
	}
	
	public SparseMatrix add(final SparseMatrix other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		final SparseMatrix ret = new SparseMatrix(this);
		for (int i = 0; i < other.sparseEntries; i++)
		{
			ret.add(other.sparseValues[i], other.sparseYs[i],
			        other.sparseXs[i]);
		}
		return ret;
	}
	
	public SparseMatrix sub(final SparseMatrix other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		final SparseMatrix ret = new SparseMatrix(this);
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
		final Map<IntCoordinates, Double> list = new HashMap<>();
		for (int i = 0; i < sparseEntries; i++)
		{
			final IntCoordinates coordinates = new IntCoordinates(sparseYs[i], sparseXs[i]);
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
	public SparseMatrix mul(final double scalar)
	{
		final SparseMatrix ret = new SparseMatrix(this);
		for (int i = 0; i < sparseEntries; i++)
		{
			ret.sparseValues[i] *= scalar;
		}
		return ret;
	}
	
	@Override
	public List<Vector> unfoldDimension(final int dimension)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (dimension > 1)
				throw new IllegalArgumentException("Matrix is two dimensional");
		final List<Vector> ret = new ArrayList<>(getShape().get(dimension));
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
		final SparseMatrix ret = new SparseMatrix(cols, rows);
		ret.sparseEntries = sparseEntries;
		ret.sparseValues = sparseValues.clone();
		ret.sparseXs = sparseYs.clone();
		ret.sparseYs = sparseXs.clone();
		return ret;
	}
	
	@Override
	public DenseVector mvMul(final Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (vector.getLength()))
				throw new IllegalArgumentException("Incompatible sizes");
		final DenseVector ret = new DenseVector(getRows());
		for (int i = 0; i < sparseEntries; i++)
			ret.add(sparseValues[i] * vector.at(sparseXs[i]), sparseYs[i]);
		return ret;
	}
	
	@Override
	public DenseVector tvMul(final Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != (vector.getLength()))
				throw new IllegalArgumentException("Incompatible sizes");
		final DenseVector ret = new DenseVector(getCols());
		for (int i = 0; i < sparseEntries; i++)
			ret.add(sparseValues[i] * vector.at(sparseYs[i]), sparseXs[i]);
		return ret;
	}
	
	@Override
	public DirectlySolvable mmMul(final Matrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (matrix.getRows()))
				throw new IllegalArgumentException("Incompatible sizes");
		if (matrix instanceof SparseMatrix)
		{
			final Map<Integer, Map<Integer, Double>> rowEntries =
				getCoordinateEntryList().entrySet()
				                        .stream()
				                        .map(e -> new Tuple2<>(e.getKey(), e.getValue()))
				                        .collect(Collectors.groupingBy(e -> e._1.get(0),
				                                                       new MapKeySelectCollector<>(key -> key.get(
					                                                       1))));
			final Map<Integer, Map<Integer, Double>> colEntries =
				((SparseMatrix) matrix).getCoordinateEntryList()
				                       .entrySet()
				                       .stream()
				                       .map(e -> new Tuple2<>(e.getKey(), e.getValue()))
				                       .collect(Collectors.groupingBy(e -> e._1.get(1),
				                                                      new MapKeySelectCollector<>(key -> key.get(
					                                                      0))));
			final SparseMatrix ret = new SparseMatrix(getRows(), matrix.getCols());
			
			Stream<Map.Entry<Integer, Map<Integer, Double>>> str = rowEntries.entrySet()
			                                                                 .stream();
			if (PerformanceArguments.getInstance().parallelizeThreads)
				str = str.parallel();
			str.forEach(rowEntry ->
			            {
				            final int rowIndex = rowEntry.getKey();
				            final Map<Integer, Double> row = rowEntry.getValue();
				            for (final Map.Entry<Integer, Map<Integer, Double>> columnEntry : colEntries.entrySet())
				            {
					            final int colIndex = columnEntry.getKey();
					            final Map<Integer, Double> col = columnEntry.getValue();
					            double contraction = 0;
					            for (final int k : row.keySet())
						            contraction += col.getOrDefault(k, 0.0);
					            ret.add(contraction, rowIndex, colIndex);
				            }
			            });
			return ret;
			
			/*final SparseMatrix m = (SparseMatrix) matrix;
			final SparseMatrix ret = new SparseMatrix(getRows(), matrix.getCols());
			IntStream stream = IntStream.range(0, sparseEntries);
			if (PerformanceArguments.getInstance().parallelizeThreads)
				stream = stream.parallel();
			stream.forEach(i ->
			               {
				               for (int j = 0; j < m.sparseEntries; j++)
				               {
					               if (sparseXs[i] == m.sparseYs[j])
						               ret.add(sparseValues[i] * m.sparseValues[j],
						                       sparseYs[i],
						                       m.sparseXs[j]);
				               }
			               });
			return ret;*/
		} else if (matrix.isSparse())
		{
			return this.mmMul(new SparseMatrix(matrix));
		} else
			return (new DenseMatrix(this)).mmMul(matrix);
	}
	
	@Override
	public SparseMatrix getLowerTriangleMatrix()
	{
		final SparseMatrix mat = new SparseMatrix(getRows(), getCols());
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
		final SparseMatrix mat = new SparseMatrix(getRows(), getCols());
		for (int i = 0; i < sparseEntries; i++)
		{
			if (sparseYs[i] < sparseXs[i])
				mat.add(sparseValues[i], sparseYs[i], sparseXs[i]);
		}
		return mat;
	}
	
	@Override
	public boolean almostEqual(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Tensors are of different size");
		}
		if (!other.isSparse())
			return other.almostEqual(this);
		final ImmutableMap<IntCoordinates, Double> myValues = getCoordinateEntryList();
		final ImmutableMap<IntCoordinates, Double> otherValues = other.getCoordinateEntryList();
		final double absmax = absMaxElement() + other.absMaxElement();
		if (coordinateEntryListsEqual(myValues, otherValues, absmax)) return false;
		if (otherValues.size() != myValues.size())
			System.out.println("matrices have different numbers of values");
		return otherValues.size() == myValues.size();
	}
	
	public static boolean coordinateEntryListsEqual(final ImmutableMap<IntCoordinates, Double> myValues,
	                                                final ImmutableMap<IntCoordinates, Double> otherValues,
	                                                final double absmax)
	{
		for (final Map.Entry<IntCoordinates, Double> entry : myValues.entrySet())
		{
			if (otherValues.containsKey(entry.getKey()))
			{
				if (!DoubleCompare.almostEqualAfterOps(entry.getValue(),
				                                       otherValues.get(entry.getKey()),
				                                       absmax * 10,
				                                       myValues.size() + otherValues.size()))
				{
					System.out.println(entry.getValue() + " != " + otherValues.get(entry.getKey()) +
						                   " difference: " + Math.abs(
						entry.getValue() - otherValues.get(entry.getKey())) + " at " +
						                   entry.getKey());
					return true;
				}
			} else
			{
				System.out.println("Value at" + entry.getKey() + " was not found in second matrix");
				return true;
			}
		}
		return false;
	}
	
	@Override
	public boolean equals(final Object obj)
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
		final SparseMatrix mat = new SparseMatrix(getRows(), getCols());
		for (int i = 0; i < sparseEntries; i++)
		{
			if (sparseYs[i] == sparseXs[i])
				mat.add(sparseValues[i], sparseYs[i], sparseXs[i]);
		}
		return mat;
	}
}
