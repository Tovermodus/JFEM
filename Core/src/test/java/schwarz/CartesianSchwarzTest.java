package schwarz;

import linalg.*;
import mixed.QkQkFunction;
import mixed.TaylorHoodSpace;
import org.junit.Test;
import tensorproduct.ContinuousTPFESpace;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.geometry.TPCell;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class CartesianSchwarzTest
{
	public static SparseMatrix createSparseMatrix(final int n)
	{
		final Random generator = new Random(14239078);
		final SparseMatrix ret = new SparseMatrix(n, n);
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			if (generator.nextDouble() < 0.4)
				ret.add((int) (10 * generator.nextDouble() - 5), c);
			if (c.get(0) == c.get(1))
				ret.add(15, c);
		}
		return ret.add(ret.transpose());
	}
	
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
			                                new AdditiveSubspaceCorrection<>(1), new DirectSolver());
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
		final IntCoordinates cells = new IntCoordinates(24, 27);
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
			                                new AdditiveSubspaceCorrection<>(1), new DirectSolver());
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
			                                new AdditiveSubspaceCorrection<>(1), new DirectSolver());
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
	
	@Test
	public void testSchwarzProjection()
	{
		final IntCoordinates cells = new IntCoordinates(2, 2);
		final ContinuousTPFESpace space = new ContinuousTPFESpace(CoordinateVector.fromValues(0, 0),
		                                                          CoordinateVector.fromValues(1, 1),
		                                                          cells);
		space.assembleCells();
		space.assembleFunctions(1);
		final SparseMatrix A = createSparseMatrix(space.getShapeFunctions()
		                                               .size());
		
		final IntCoordinates partitions = new IntCoordinates(2, 2);
		final int overlap = 0;
		final CartesianUpFrontSchwarz<ContinuousTPShapeFunction> schwarz
			= new CartesianUpFrontSchwarz<>(A,
			                                space,
			                                partitions,
			                                overlap,
			                                new AdditiveSubspaceCorrection<>(1), new DirectSolver());
		for (int i = 0; i < schwarz.getPatchCount(); i++)
		{
			final Matrix R = schwarz.getRestrictionOperator(i);
			final Matrix Aiinv = new DenseMatrix(schwarz.getLocalOperator(i)).inverse();
			final Matrix P
				= R.tmMul(Aiinv)
				   .mmMul(R)
				   .mmMul(A);
			assertTrue(P.mmMul(P)
			            .sub(P)
			            .absMaxElement() < 1e-14);
			assertTrue(P.tmMul(A)
			            .sub(A.mmMul(P))
			            .absMaxElement() < 1e-14);
		}
	}
	
	@Test
	public void testRestrictionMatrix()
	{
		RestrictionMatrix r = new RestrictionMatrix(5, 10);
		r.addSmallMatrixInPlaceAt(SparseMatrix.identity(5), 0, 0);
		final SparseMatrix A = createSparseMatrix(10);
//		System.out.println(A);
//		System.out.println(r.selectFrom(A));
//		System.out.println(r.mmMul(A)
//		                    .mtMul(r));
		assertEquals(r.selectFrom(A),
		             r.mmMul(A)
		              .mtMul(r));
		r = new RestrictionMatrix(5, 10);
		r.add(1, 0, 0);
		r.add(1, 1, 2);
		r.add(1, 2, 1);
		r.add(1, 3, 9);
		r.add(1, 4, 8);
//		System.out.println(A);
//		System.out.println(r.selectFrom(A));
//		System.out.println(r.mmMul(A)
//		                    .mtMul(r));
		assertEquals(r.selectFrom(A),
		             r.mmMul(A)
		              .mtMul(r));
		r = new RestrictionMatrix(5, 10);
		r.add(1, 1, 0);
		r.add(1, 3, 2);
		r.add(1, 2, 1);
		r.add(1, 4, 9);
		r.add(1, 0, 8);
//		System.out.println(A);
//		System.out.println(r.selectFrom(A));
//		System.out.println(r.mmMul(A)
//		                    .mtMul(r));
		assertEquals(r.selectFrom(A),
		             r.mmMul(A)
		              .mtMul(r));
		r = new RestrictionMatrix(5, 10);
		r.add(1, 0, 4);
		r.add(1, 1, 3);
		r.add(1, 2, 2);
		r.add(1, 3, 5);
		r.add(1, 4, 8);
//		System.out.println(r);
//		System.out.println(A);
//		System.out.println(r.selectFrom(A));
//		System.out.println(r.mmMul(A)
//		                    .mtMul(r));
		assertEquals(r.selectFrom(A),
		             r.mmMul(A)
		              .mtMul(r));
	}
}
