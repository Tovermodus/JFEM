package linalg;

import basic.PerformanceArguments;

import java.util.ArrayList;
import java.util.Map;
import java.util.stream.IntStream;

public class SparseMvMul
	implements VectorMultiplyable
{
	private final int rows;
	private final int cols;
	final ArrayList<double[]> sparseRowValues;
	final ArrayList<double[]> sparseColValues;
	final ArrayList<int[]> sparseRowXs;
	final ArrayList<int[]> sparseColYs;
	boolean parallel;
	
	public SparseMvMul(final Matrix s)
	{
		parallel = PerformanceArguments.getInstance().parallelizeThreads;
		rows = s.getRows();
		cols = s.getCols();
		sparseRowValues = new ArrayList<>();
		sparseColValues = new ArrayList<>();
		sparseRowXs = new ArrayList<>();
		sparseColYs = new ArrayList<>();
		final Map<IntCoordinates, Double> coordinateEntries = s.getCoordinateEntryList();
		final int[] rowCounts = new int[rows];
		final int[] colCounts = new int[cols];
		for (final IntCoordinates i : coordinateEntries.keySet())
		{
			rowCounts[i.get(0)]++;
			colCounts[i.get(1)]++;
		}
		for (int i = 0; i < rows; i++)
		{
			sparseRowValues.add(new double[rowCounts[i]]);
			sparseRowXs.add(new int[rowCounts[i]]);
		}
		for (int i = 0; i < cols; i++)
		{
			sparseColValues.add(new double[colCounts[i]]);
			sparseColYs.add(new int[colCounts[i]]);
		}
		s.getCoordinateEntryList()
		 .forEach((k, v) ->
		          {
			          final int Y = k.get(0);
			          final int X = k.get(1);
			          rowCounts[Y]--;
			          colCounts[X]--;
			          sparseColYs.get(X)[colCounts[X]] = Y;
			          sparseColValues.get(X)[colCounts[X]] = v;
			          sparseRowXs.get(Y)[rowCounts[Y]] = X;
			          sparseRowValues.get(Y)[rowCounts[Y]] = v;
		          });
	}
	
	@Override
	public int getVectorSize()
	{
		return cols;
	}
	
	@Override
	public int getTVectorSize()
	{
		return rows;
	}
	
	@Override
	public DenseVector mvMul(final Vector vector)
	{
		final double[] vectorVals = vector.asArray();
		final double[] ret = new double[rows];
		IntStream str = IntStream.range(0, rows);
		if (parallel)
			str = str.parallel();
		str.forEach(row ->
		            {
			            final double[] sparseVals = sparseRowValues.get(row);
			            final int[] sparseXs = sparseRowXs.get(row);
			            double ret_component = 0;
			            for (int i = 0; i < sparseVals.length; i++)
				            ret_component += sparseVals[i] * vectorVals[sparseXs[i]];
			            ret[row] = ret_component;
		            });
		return new DenseVector(ret, true);
	}
	
	@Override
	public DenseVector tvMul(final Vector vector)
	{
		final double[] vectorVals = vector.asArray();
		final double[] ret = new double[cols];
		IntStream str = IntStream.range(0, cols);
		if (parallel)
			str = str.parallel();
		str.forEach(col ->
		            {
			            final double[] sparseVals = sparseColValues.get(col);
			            final int[] sparseYs = sparseColYs.get(col);
			            double ret_component = 0;
			            for (int i = 0; i < sparseVals.length; i++)
				            ret_component += sparseVals[i] * vectorVals[sparseYs[i]];
			            ret[col] = ret_component;
		            });
		return new DenseVector(ret, true);
	}
}
