package schwarz;

import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import linalg.IntCoordinates;
import linalg.Matrix;
import linalg.SparseMatrix;

public class RestrictionMatrix
	extends SparseMatrix
{
	public RestrictionMatrix(final int rows, final int cols, final int size)
	{
		super(rows, cols, size);
	}
	
	public RestrictionMatrix(final int rows, final int cols)
	{
		super(rows, cols);
	}
	
	public SparseMatrix selectFrom(final SparseMatrix globalMatrix)
	{
		final SparseMatrix ret = new SparseMatrix(getRows(), getRows());
		final Int2IntMap selectedIndices = new Int2IntOpenHashMap();
		selectedIndices.defaultReturnValue(-1);
		for (int i = 0; i < sparseEntries; i++)
			if (sparseValues[i] != 0)
				selectedIndices.put(sparseXs[i], sparseYs[i]);
		globalMatrix.getCoordinateEntryList()
		            .entrySet()
		            .stream()
		            .parallel()
		            .forEach(entry ->
		                     {
			                     final IntCoordinates coord = entry.getKey();
			                     final int newY = selectedIndices.get(coord.get(0));
			                     final int newX = selectedIndices.get(coord.get(1));
			                     if (newX != -1 && newY != -1)
				                     ret.add(entry.getValue(), newY, newX);
		                     });
		return ret;
	}
	
	public SparseMatrix selectFrom(final Matrix globalMatrix)
	{
		return selectFrom(new SparseMatrix(globalMatrix));
	}
}
