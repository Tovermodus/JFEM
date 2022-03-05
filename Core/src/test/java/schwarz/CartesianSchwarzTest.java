package schwarz;

import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.IntCoordinates;
import linalg.SparseMatrix;
import mixed.QkQkFunction;
import mixed.TaylorHoodSpace;
import org.junit.Test;
import tensorproduct.geometry.TPCell;

import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class CartesianSchwarzTest
{
	@Test
	public void testPatches()
	{
		final IntCoordinates cells = new IntCoordinates(4, 6);
		final TaylorHoodSpace space = new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
		                                                  CoordinateVector.fromValues(1, 1),
		                                                  cells);
		space.assembleCells();
		space.assembleFunctions(1);
		final SparseMatrix s = new SparseMatrix(space.getShapeFunctions()
		                                             .size(),
		                                        space.getShapeFunctions()
		                                             .size());
		final IntCoordinates partitions = new IntCoordinates(2, 3);
		final int overlap = 0;
		final CartesianUpFrontSchwarz<QkQkFunction> schwarz
			= new CartesianUpFrontSchwarz<>(s,
			                                space,
			                                partitions,
			                                overlap,
			                                new AdditiveSubspaceCorrection<>(1));
		assertEquals(schwarz.cellPatches.size(), partitions.size());
		assertEquals(schwarz.getPatchCount(), partitions.size());
		final Set<TPCell> allCells = new HashSet<>();
		for (int i = 0; i < schwarz.getPatchCount(); i++)
		{
			assertEquals((long) cells.get(0) / partitions.get(0) * cells.get(1) / partitions.get(1),
			             schwarz.getCellPatch(i)
			                    .size());
			for (final TPCell c : schwarz.getCellPatch(i))
				assertTrue(allCells.add(c));
		}
		assertEquals(allCells.size(), cells.size());
		space.getCells()
		     .forEach(allCells::remove);
		assertEquals(0, allCells.size());
	}
	
	@Test
	public void testPatches2()
	{
		final IntCoordinates cells = new IntCoordinates(44, 99);
		final TaylorHoodSpace space = new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
		                                                  CoordinateVector.fromValues(1, 1),
		                                                  cells);
		space.assembleCells();
		space.assembleFunctions(1);
		final SparseMatrix s = new SparseMatrix(space.getShapeFunctions()
		                                             .size(),
		                                        space.getShapeFunctions()
		                                             .size());
		final IntCoordinates partitions = new IntCoordinates(4, 9);
		final int overlap = 0;
		final CartesianUpFrontSchwarz<QkQkFunction> schwarz
			= new CartesianUpFrontSchwarz<>(s,
			                                space,
			                                partitions,
			                                overlap,
			                                new AdditiveSubspaceCorrection<>(1));
		assertEquals(schwarz.cellPatches.size(), partitions.size());
		assertEquals(schwarz.getPatchCount(), partitions.size());
		final Set<TPCell> allCells = new HashSet<>();
		for (int i = 0; i < schwarz.getPatchCount(); i++)
		{
			assertEquals((long) cells.get(0) / partitions.get(0) * cells.get(1) / partitions.get(1),
			             schwarz.getCellPatch(i)
			                    .size());
			for (final TPCell c : schwarz.getCellPatch(i))
				assertTrue(allCells.add(c));
		}
		assertEquals(allCells.size(), cells.size());
		space.getCells()
		     .forEach(allCells::remove);
		assertEquals(0, allCells.size());
	}
	
	@Test
	public void testRestriction()
	{
		final IntCoordinates cells = new IntCoordinates(6, 6);
		final TaylorHoodSpace space = new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
		                                                  CoordinateVector.fromValues(1, 1),
		                                                  cells);
		space.assembleCells();
		space.assembleFunctions(1);
		final SparseMatrix s = new SparseMatrix(space.getShapeFunctions()
		                                             .size(),
		                                        space.getShapeFunctions()
		                                             .size());
		final IntCoordinates partitions = new IntCoordinates(2, 3);
		final int overlap = 0;
		final CartesianUpFrontSchwarz<QkQkFunction> schwarz
			= new CartesianUpFrontSchwarz<>(s,
			                                space,
			                                partitions,
			                                overlap,
			                                new AdditiveSubspaceCorrection<>(1));
		DenseVector v = new DenseVector(space.getShapeFunctions()
		                                     .size());
		for (int i = 0; i < space.getShapeFunctions()
		                         .size(); i++)
			v.set(i, i);
		for (int i = 0; i < schwarz.getPatchCount(); i++)
		{
			v = v.sub(schwarz.getGlobalVector(i, schwarz.getLocalVector(i, v)));
		}
		assertTrue(v.euclidianNorm() < 1e-14);
	}
}
