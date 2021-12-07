package distorted.geometry;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableCollection;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import io.vavr.Tuple2;
import kotlin.Pair;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.*;

public interface DistortedGrid
{
	default Pair<ImmutableSet<DistortedCell>, ImmutableSet<DistortedFace>> createGrid(final int refinements)
	{
		
		Multimap<DistortedCell, DistortedCell> refinedCells = generateCells();
		Collection<DistortedFace> genFaces = generateFaces(refinedCells);
		for (int i = 0; i < refinements; i++)
		{
			System.out.println("refine");
			refinedCells = refineCells(refinedCells);
			genFaces = refineFaces(refinedCells);
			for (final DistortedFace f : genFaces)
			{
				if (!f.isBoundaryFace() && f.getCells()
				                            .size() != 2)
					System.out.println(f);
			}
		}
		int hash = 379;
		final ImmutableSet<DistortedCell> cells = ImmutableSet.copyOf(refinedCells.values());
		System.out.println(cells.size());
		final ImmutableSet<DistortedFace> faces = ImmutableSet.copyOf(genFaces);
		for (final DistortedCell cell : cells)
			cell.setDone(hash++);
		return new Pair<>(cells, faces);
	}
	
	static HashSet<DistortedCell> getCellsBelongingToFace(final DistortedFace f,
	                                                      final Collection<DistortedCell> from)
	{
		final HashSet<DistortedCell> cellsBelonging = new HashSet<>();
		for (final DistortedCell cell : from)
		{
			if (f.isOnCell(cell))
				cellsBelonging.add(cell);
		}
		return cellsBelonging;
	}
	
	Multimap<DistortedCell, DistortedCell> generateCells();
	
	default Collection<DistortedFace> generateFaces(final Multimap<DistortedCell, DistortedCell> genCells)
	{
		final HashSet<DistortedFace> faces = new HashSet<>();
		for (final DistortedCell cell : genCells.values())
		{
			createFaces(faces, cell, genCells.values());
		}
		return faces;
	}
	
	default HashSet<DistortedFace> refineFaces(final Multimap<DistortedCell, DistortedCell> refinedCells)
	{
		final HashSet<DistortedFace> faces = new HashSet<>();
		for (final DistortedCell cell : refinedCells.keySet())
		{
			final List<DistortedCell> refinedNeighbouringCells = getRefinedNeighbouringCells(refinedCells,
			                                                                                 cell);
			for (final DistortedCell refinedCell : refinedCells.get(cell))
			{
				createFaces(faces, refinedCell, refinedNeighbouringCells);
			}
		}
		return faces;
	}
	
	void createFaces(final HashSet<DistortedFace> faces, final DistortedCell cell,
	                 final Collection<DistortedCell> neighbouringCells);
	
	@NotNull
	default List<DistortedCell> getRefinedNeighbouringCells(final Multimap<DistortedCell, DistortedCell> refinedCells,
	                                                        final DistortedCell cell)
	{
		final List<DistortedCell> neighbouringCells = new ArrayList<>();
		for (final DistortedFace f : cell.faces)
			neighbouringCells.addAll(f.getCells());
		final List<DistortedCell> refinedNeighbouringCells = new ArrayList<>();
		for (final DistortedCell c : neighbouringCells)
			refinedNeighbouringCells.addAll(refinedCells.get(c));
		return refinedNeighbouringCells;
	}
	
	default CoordinateVector mean(final CoordinateVector... points)
	{
		return Arrays.stream(points)
		             .reduce(points[0].mul(0), CoordinateVector::add)
		             .mul(1. / points.length);
	}
	
	List<DistortedCell> partitionCell(final DistortedCell cell);
	
	default Multimap<DistortedCell, DistortedCell> refineCells(final Multimap<DistortedCell, DistortedCell> genCells)
	{
		final Multimap<DistortedCell, DistortedCell> refinedCells = HashMultimap.create();
		for (final DistortedCell cell : genCells.values())
			refinedCells.putAll(cell, partitionCell(cell));
		return refinedCells;
	}
	
	List<CoordinateVector> generatePlotPoints(final int resolution);
	
	List<CoordinateVector> generateIsotropicPlotPoints(final int resolution);
	
	List<Tuple2<DistortedCell, CoordinateVector>> generateReferencePlotPoints(final int resolution);
	
	ImmutableCollection<DistortedFace> getFaces();
	
	ImmutableCollection<DistortedCell> getCells();
	
	int getDimension();
}
