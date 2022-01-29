package linalg;

import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.*;

public class TriangleMatrixTest
{
	public static LowerTriangularMatrix getLowerTriangle()
	{
		final Random generator = new Random(9831415);
		final DenseMatrix dm = new DenseMatrix(10, 10);
		for (final IntCoordinates c : dm.getShape()
		                                .range())
			dm.add(generator.nextDouble(), c);
		return new LowerTriangularMatrix(dm);
	}
	
	public static UpperTriangularMatrix getUpperTriangle()
	{
		
		final Random generator = new Random(98314315);
		final DenseMatrix dm = new DenseMatrix(10, 10);
		for (final IntCoordinates c : dm.getShape()
		                                .range())
			dm.add(generator.nextDouble(), c);
		return new UpperTriangularMatrix(dm);
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
		final LowerTriangularMatrix l = getLowerTriangle();
		final UpperTriangularMatrix u = getUpperTriangle();
		assertThrows(IllegalArgumentException.class, () -> l.add(1, 1, 3));
		assertThrows(IllegalArgumentException.class, () -> u.add(1, 3, 1));
	}
	
	@Test
	public void testInsert()
	{
		final LowerTriangularMatrix l = getLowerTriangle();
		final UpperTriangularMatrix u = getUpperTriangle();
		final DenseMatrix dl = new DenseMatrix(l);
		final DenseMatrix du = new DenseMatrix(u);
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
		assertTrue(getLowerTriangle().transpose() instanceof UpperTriangularMatrix);
		assertEquals(getLowerTriangle().transpose()
		                               .transpose(), getLowerTriangle());
		assertTrue(getUpperTriangle().transpose() instanceof LowerTriangularMatrix);
		assertEquals(getUpperTriangle().transpose()
		                               .transpose(), getUpperTriangle());
	}
	
	@Test
	public void testSolve()
	{
		final LowerTriangularMatrix l = getLowerTriangle();
		assertEquals(l.solve(getRhs(l)), new DenseMatrix(l).solve(getRhs(l)));
		final UpperTriangularMatrix u = getUpperTriangle();
		assertEquals(u.solve(getRhs(u)), new DenseMatrix(u).solve(getRhs(u)));
	}
	
	@Test
	public void testInverse()
	{
		final LowerTriangularMatrix l = getLowerTriangle();
		assertEquals(l.inverse(), new DenseMatrix(l).inverse());
		final UpperTriangularMatrix u = getUpperTriangle();
		assertEquals(u.inverse(), new DenseMatrix(u).inverse());
	}
}
