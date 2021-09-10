package linalg;

import com.google.common.collect.ImmutableMap;
import org.junit.Test;

import java.util.stream.IntStream;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class DenseVectorTest
{
	private final int n = 10000;
	
	private DenseVector createVector()
	{
		final double[] values = new double[n];
		for (int i = 0; i < n; i++)
		{
			values[i] = i;
		}
		return new DenseVector(values);
	}
	
	@Test
	public void testGetters()
	{
		final DenseVector vector = createVector();
		final DenseVector vector2 = createVector();
		assertEquals(vector.getLength(), n);
		assertEquals(vector.getOrder(), 1);
		assertEquals(vector.getShape(), new IntCoordinates(n));
		final ImmutableMap<IntCoordinates, Double> entries = vector.getCoordinateEntryList();
		for (int i = 1; i < n; i++)
		{
			assertEquals(entries.get(new IntCoordinates(i)).doubleValue(), i);
		}
		assertEquals(vector.size(), n);
		assertEquals(vector, vector2);
	}
	
	@Test
	public void testAtSetAdd()
	{
		final DenseVector vector = createVector();
		for (int i = 0; i < n; i++)
			assertEquals(vector.at(i), i);
		IntStream.range(0, n).parallel().forEach(i -> vector.add(7.7, i));
		for (int i = 0; i < n; i++)
			assertEquals(vector.at(i), 7.7 + i);
		assertTrue(IntStream.range(0, n).parallel().allMatch(i -> vector.at(i) == 7.7 + i));
		IntStream.range(0, n).parallel().forEach(i -> vector.add(-7.7, i));
		assertTrue(IntStream.range(0, n).parallel().allMatch(i -> Math.abs(vector.at(i) - i) < n * 1e-15));
	}
	
	@Test
	public void testAddInPlace()
	{
		final DenseVector vector = createVector();
		for (int i = 0; i < n; i++)
			assertEquals(vector.at(i), i);
		final DenseVector v2 = createVector();
		v2.addInPlace(vector);
		IntStream.range(0, n).parallel().forEach(i -> vector.set(vector.at(i) * 2, i));
		assertTrue(vector.almostEqual(v2));
		assertEquals(vector, v2);
	}
	
	@Test
	public void testMulInPlace()
	{
		final DenseVector vector = createVector();
		for (int i = 0; i < n; i++)
			assertEquals(vector.at(i), i);
		final DenseVector v2 = createVector();
		v2.mulInPlace(17);
		IntStream.range(0, n).parallel().forEach(i -> vector.set(vector.at(i) * 17, i));
		assertTrue(vector.almostEqual(v2));
		assertEquals(vector, v2);
	}
	
	@Test
	public void testAddVector()
	{
		final DenseVector vector = createVector();
		for (int i = 0; i < n; i++)
			assertEquals(vector.at(i), i);
		final DenseVector v2 = createVector();
		IntStream.range(0, n).parallel().forEach(i -> vector.set(vector.at(i) * 2, i));
		assertTrue(vector.almostEqual(v2.add(v2)));
		assertEquals(vector, v2.add(v2));
	}
	
	@Test
	public void testMul()
	{
		final DenseVector vector = createVector();
		for (int i = 0; i < n; i++)
			assertEquals(vector.at(i), i);
		final DenseVector v2 = vector.mul(17);
		IntStream.range(0, n).parallel().forEach(i -> vector.set(vector.at(i) * 17, i));
		assertTrue(vector.almostEqual(v2));
		assertEquals(vector, v2);
	}
	
	@Test
	public void testOuter()
	{
		final DenseVector xv = DenseVector.vectorFromValues(2, 3, 2);
		final DenseVector yv = DenseVector.vectorFromValues(9, 8, 9);
		final DenseMatrix comp = DenseMatrix.squareMatrixFromValues(2 * 9, 3 * 9, 2 * 9,
		                                                            2 * 8, 3 * 8, 2 * 8, 2 * 9, 3 * 9, 2 * 9);
		assertEquals(comp, yv.outer(xv));
		assertTrue(comp.almostEqual(yv.outer(xv)));
		assertTrue(comp.almostEqual(xv.outer(yv).transpose()));
	}
	
	@Test
	public void testInner()
	{
		final DenseVector vector = createVector();
		final double val = vector.inner(vector.mul(0.0002));
		final double valcomp = IntStream.range(0, n).mapToDouble((i) -> i * i * 0.0002).sum();
		assertTrue(Math.abs(val - valcomp) < 1e-15 * n * n);
	}
	
	@Test
	public void testMatrixVectorProduct()
	{
		final DenseVector v = createVector();
		final SparseMatrix s = SparseMatrix.identity(n);
		assertEquals(s.mvMul(v), v);
		s.mulInPlace(10);
		assertEquals(s.mvMul(v), v.mul(10));
	}
}
