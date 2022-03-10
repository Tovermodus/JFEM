package schwarz;

import basic.ShapeFunction;
import it.unimi.dsi.fastutil.ints.IntRBTreeSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import linalg.IntCoordinates;
import linalg.Matrix;
import mixed.TaylorHoodSpace;
import tensorproduct.CartesianGridSpace;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;

public class ColoredCartesianSchwarz<ST extends ShapeFunction<TPCell, TPFace, ?, ?, ?>>
	extends UpFrontSchwarz<TPCell, TPFace>
{
	final List<Set<TPCell>> cellPatches;
	private final CartesianGridSpace<ST, ?, ?, ?> space;
	final IntCoordinates partitions;
	final int overlap;
	
	public ColoredCartesianSchwarz(final Matrix globalMatrix,
	                               final CartesianGridSpace<ST, ?, ?, ?> space,
	                               final IntCoordinates partitions,
	                               final int overlap,
	                               final SystemSolver<Matrix> solver)
	{
		super(globalMatrix, new ColoredMultiplicativeSubspaceCorrection<>((TaylorHoodSpace) space), solver);
		this.space = space;
		this.partitions = partitions;
		this.overlap = overlap;
		if (space.getDimension() != 2)
			throw new IllegalArgumentException("Only in 2D");
		final IntCoordinates underlyingCells = space.grid.cellsPerDimension;
		if (underlyingCells.get(0) % partitions.get(0) != 0)
			throw new IllegalArgumentException("First dimension does not fit");
		if (underlyingCells.get(1) % partitions.get(1) != 0)
			throw new IllegalArgumentException("Second dimension does not fit");
		final IntCoordinates cellsPerPartition = new IntCoordinates(underlyingCells.get(0) / partitions.get(0),
		                                                            underlyingCells.get(1) / partitions.get(1));
		System.out.println("Partitions " + partitions + " Cells " + underlyingCells + " Cells per partition " + cellsPerPartition + " Overlap " + overlap);
		if (overlap >= cellsPerPartition.get(0) || overlap >= cellsPerPartition.get(1))
			throw new IllegalArgumentException("Too Much Overlap, can't color. Max overlap is " + Math.min(
				cellsPerPartition.get(0) - 1,
				cellsPerPartition.get(1) - 1));
		cellPatches = new ArrayList<>();
		for (final IntCoordinates partition : partitions.range())
		{
			final HashSet<TPCell> patch = new HashSet<>();
			for (int i = -overlap; i < cellsPerPartition.get(0) + overlap; i++)
			{
				for (int j = -overlap; j < cellsPerPartition.get(1) + overlap; j++)
				{
					final IntCoordinates cellCoords =
						new IntCoordinates(partition.get(0) * cellsPerPartition.get(0) + i,
						                   partition.get(1) * cellsPerPartition.get(1) + j);
					if (cellCoords.get(0) >= 0 && cellCoords.get(1) >= 0
						&& cellCoords.compareTo(underlyingCells) < 0)
					{
						patch.add(space.grid.cellsByCoordinates.get(cellCoords));
					}
				}
			}
			cellPatches.add(patch);
		}
		System.out.println("partitioned");
		((ColoredMultiplicativeSubspaceCorrection<Matrix>) getSubspaceCorrection()).setColors(getColors());
		build();
	}
	
	private List<IntSet> getColors()
	{
		final List<IntSet> colors = new ArrayList<>();
		for (int i = 0; i < 9; i++)
			colors.add(new IntRBTreeSet());
		int partitionIndex = 0;
		for (final IntCoordinates partition : partitions.range())
		{
			if (partition.get(0) % 3 == 0)
			{
				if (partition.get(1) % 3 == 0)
					colors.get(0)
					      .add(partitionIndex);
				else if (partition.get(1) % 3 == 1)
					colors.get(1)
					      .add(partitionIndex);
				else
					colors.get(2)
					      .add(partitionIndex);
			} else if (partition.get(0) % 3 == 1)
			{
				if (partition.get(1) % 3 == 0)
					colors.get(3)
					      .add(partitionIndex);
				else if (partition.get(1) % 3 == 1)
					colors.get(4)
					      .add(partitionIndex);
				else
					colors.get(5)
					      .add(partitionIndex);
			} else
			{
				if (partition.get(1) % 3 == 0)
					colors.get(6)
					      .add(partitionIndex);
				else if (partition.get(1) % 3 == 1)
					colors.get(7)
					      .add(partitionIndex);
				else
					colors.get(8)
					      .add(partitionIndex);
			}
			partitionIndex++;
		}
		return colors;
	}
	
	public Collection<TPCell> getCellPatch(final int patch)
	{
		return cellPatches.get(patch);
	}
	
	@Override
	public int getPatchCount()
	{
		return cellPatches.size();
	}
	
	@Override
	public RestrictionMatrix buildRestrictionMatrix(final int patch)
	{
		final Collection<TPCell> patchCells = getCellPatch(patch);
		final Set<ST> functions = new TreeSet<ST>(Comparator.comparingInt(st -> st.getGlobalIndex()));
		for (final TPCell c : patchCells)
			functions.addAll(space.getCellSupportMapping()
			                      .get(c));
		final RestrictionMatrix s = new RestrictionMatrix(functions.size(),
		                                                  space.getShapeFunctions()
		                                                       .size());
		int i = 0;
		for (final ST f : functions)
		{
			s.add(1, i++, f.getGlobalIndex());
		}
		return s;
	}
}
