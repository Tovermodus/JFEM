package linalg;

import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.*;

public class SparseTriangleMatrixTest
{
	public static LowerTriangularSparseMatrix getLowerTriangle()
	{
		final Random generator = new Random(9831415);
		final SparseMatrix dm = new SparseMatrix(10, 10);
		for (final IntCoordinates c : dm.getShape()
		                                .range())
			if (generator.nextDouble() < 0.2)
				dm.add(generator.nextDouble(), c);
		for (int i = 0; i < 10; i++)
			dm.add(1, i, i);
		return new LowerTriangularSparseMatrix(dm);
	}
	
	public static UpperTriangularSparseMatrix getUpperTriangle()
	{
		
		final Random generator = new Random(98314315);
		final SparseMatrix dm = new SparseMatrix(10, 10);
		for (final IntCoordinates c : dm.getShape()
		                                .range())
			if (generator.nextDouble() < 0.2)
				dm.add(generator.nextDouble(), c);
		for (int i = 0; i < 10; i++)
			dm.add(1, i, i);
		return new UpperTriangularSparseMatrix(dm);
	}
	
	public static DenseVector getRhs(final Matrix m)
	{
		
		final Random generator = new Random(94315);
		final DenseVector d = new DenseVector(m.getVectorSize());
		for (final IntCoordinates c : d.getShape()
		                               .range())
			d.add(generator.nextDouble(), c);
		return d;
	}
	
	@Test
	public void testWrongInsert()
	{
		final LowerTriangularSparseMatrix l = getLowerTriangle();
		final UpperTriangularSparseMatrix u = getUpperTriangle();
		assertThrows(IllegalArgumentException.class, () -> l.add(1, 1, 3));
		assertThrows(IllegalArgumentException.class, () -> u.add(1, 3, 1));
	}
	
	@Test
	public void testInsert()
	{
		final LowerTriangularSparseMatrix l = getLowerTriangle();
		final UpperTriangularSparseMatrix u = getUpperTriangle();
		final SparseMatrix dl = new SparseMatrix(l);
		final SparseMatrix du = new SparseMatrix(u);
		l.add(1, 0, 0);
		dl.add(1, 0, 0);
		u.add(-3, 7, 7);
		du.add(-3, 7, 7);
		l.add(.1, 6, 0);
		dl.add(.1, 6, 0);
		u.add(-3.23, 7, 9);
		du.add(-3.23, 7, 9);
		assertEquals(u, du);
		assertEquals(l, dl);
	}
	
	@Test
	public void testTranspose()
	{
		assertTrue(getLowerTriangle().transpose() instanceof UpperTriangularSparseMatrix);
		assertEquals(getLowerTriangle().transpose()
		                               .transpose(), getLowerTriangle());
		assertTrue(getUpperTriangle().transpose() instanceof LowerTriangularSparseMatrix);
		assertEquals(getUpperTriangle().transpose()
		                               .transpose(), getUpperTriangle());
	}
	
	@Test
	public void testSolve()
	{
		final LowerTriangularSparseMatrix l = getLowerTriangle();
		assertEquals(l.solve(getRhs(l)), new SparseMatrix(l).solve(getRhs(l)));
		final UpperTriangularSparseMatrix u = getUpperTriangle();
		assertEquals(u.solve(getRhs(u)), new SparseMatrix(u).solve(getRhs(u)));
	}
	
	@Test
	public void testInverse()
	{
		final LowerTriangularSparseMatrix l = getLowerTriangle();
		assertEquals(l.inverse(), new SparseMatrix(l).inverse());
		final UpperTriangularSparseMatrix u = getUpperTriangle();
		assertEquals(u.inverse(), new SparseMatrix(u).inverse());
	}
}
