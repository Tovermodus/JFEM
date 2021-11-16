package linalg;

import basic.PerformanceArguments;

public interface MutableMatrix
	extends Matrix, MutableTensor
{
	void deleteRow(int row);
	
	void deleteColumn(int column);
	
	default void addSmallMatrixInPlaceAt(final Matrix small, final int... coordinates)
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
		for (final IntCoordinates c : small.getShape()
		                                   .range())
		{
			add(small.at(c), c.get(0) + coordinates[0], c.get(1) + coordinates[1]);
		}
	}
	
	default void addColumn(final Vector vector, final int column)
	{
		addColumn(vector, column, 0);
	}
	
	default void addColumn(final Vector vector, final int column, final int start)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (start + vector.getLength() > getRows())
				throw new IllegalArgumentException("small Vector too large for position");
		}
		for (final IntCoordinates c : vector.getShape()
		                                    .range())
		{
			if (vector.at(c) != 0)
				add(vector.at(c), c.get(0) + start, column);
		}
	}
	
	default void addRow(final Vector vector, final int row)
	{
		addRow(vector, row, 0);
	}
	
	default void addRow(final Vector vector, final int row, final int start)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (start + vector.getLength() > getCols())
				throw new IllegalArgumentException("small Vector too large for position");
		}
		for (final IntCoordinates c : vector.getShape()
		                                    .range())
		{
			add(vector.at(c), row, c.get(0) + start);
		}
	}
}
