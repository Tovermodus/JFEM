package schwarz;

import basic.ShapeFunction;
import linalg.DenseVector;
import linalg.IntCoordinates;
import linalg.Matrix;
import linalg.Vector;
import tensorproduct.CartesianGridSpace;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;
import java.util.stream.IntStream;

public class CartesianUpFrontSchwarz<ST extends ShapeFunction<TPCell, TPFace, ?, ?, ?>>
	extends UpFrontSchwarz<TPCell, TPFace>
{
	final List<Set<TPCell>> cellPatches;
	final SubspaceCorrection<Matrix> subspaceCorrection;
	private final CartesianGridSpace<ST, ?, ?, ?> space;
	
	public CartesianUpFrontSchwarz(final Matrix globalMatrix,
	                               final CartesianGridSpace<ST, ?, ?, ?> space,
	                               final IntCoordinates partitions,
	                               final int overlap,
	                               final SubspaceCorrection<Matrix> subspaceCorrection,
	                               final SystemSolver<Matrix> solver)
	{
		super(globalMatrix, subspaceCorrection, solver);
		this.space = space;
		if (space.getDimension() != 2)
			throw new IllegalArgumentException("Only in 2D");
		final IntCoordinates underlyingCells = space.grid.cellsPerDimension;
		final IntCoordinates cellsPerPartition;
		this.subspaceCorrection = subspaceCorrection;
		cellPatches = new ArrayList<>();
		if (partitions.get(0) != -1)
		{
			if (underlyingCells.get(0) % partitions.get(0) != 0)
				throw new IllegalArgumentException("First dimension does not fit");
			if (underlyingCells.get(1) % partitions.get(1) != 0)
				throw new IllegalArgumentException("Second dimension does not fit");
			cellsPerPartition = new IntCoordinates(underlyingCells.get(0) / partitions.get(0),
			                                       underlyingCells.get(1) / partitions.get(1));
			
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
							&& cellCoords.get(0) < underlyingCells.get(0)
							&& cellCoords.get(1) < underlyingCells.get(1))
						{
							final TPCell cell
								= space.grid.cellsByCoordinates.get(cellCoords);
							patch.add(cell);
						}
					}
				}
				cellPatches.add(patch);
			}
		} else
		{
			final IntCoordinates offset = new IntCoordinates(1, 1);
			for (final IntCoordinates partition : underlyingCells.sub(offset)
			                                                     .range())
			{
				final IntCoordinates upperleftCorner = offset.add(partition);
				System.out.println("p " + upperleftCorner);
				final HashSet<TPCell> patch = new HashSet<>();
				for (int i = -overlap; i < overlap; i++)
				{
					for (int j = -overlap; j < overlap; j++)
					{
						final IntCoordinates cellCoords =
							new IntCoordinates(upperleftCorner.get(0) + i,
							                   upperleftCorner.get(1) + j);
						if (cellCoords.get(0) >= 0 && cellCoords.get(1) >= 0
							&& cellCoords.get(0) < underlyingCells.get(0)
							&& cellCoords.get(1) < underlyingCells.get(1))
						{
							final TPCell cell
								= space.grid.cellsByCoordinates.get(cellCoords);
							patch.add(cell);
						}
					}
				}
				cellPatches.add(patch);
			}
		}
		System.out.println("partitioned");
		build();
		
		IntStream.range(0, getPatchCount())
		         .forEach(i ->
		                  {
			                  final Vector v = new DenseVector(globalMatrix.getVectorSize());
			                  final linalg.Vector localRes = getLocalVector(i,
			                                                                v);
			                  final Vector localSol = solveLocalSystem(i,
			                                                           localRes);
		                  });
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
	public SubspaceCorrection<Matrix> getSubspaceCorrection()
	{
		return subspaceCorrection;
	}
	
	@Override
	public RestrictionMatrix buildRestrictionMatrix(final int patch)
	{
		final Collection<TPCell> patchCells = getCellPatch(patch);
		final Set<ST> functions = new TreeSet<ST>(Comparator.comparingInt(st -> st.getGlobalIndex()));
		for (final TPCell c : patchCells)
		{
			functions.addAll(space.getShapeFunctionsWithSupportOnCell(c));
		}
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
