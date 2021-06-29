package linalg;

import com.google.common.collect.ImmutableMap;
import  org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class ConcurrentDenseVectorTest
{
	private final int n = 100;
	private ConcurrentDenseVector createVector()
	{
		double[] values = new double[n];
		for (int i = 0; i < n; i++)
		{
			values[i] = i;
		}
		return new ConcurrentDenseVector(values);
	}
	@Test
	public void testGetters() {
		ConcurrentDenseVector vector = createVector();
		ConcurrentDenseVector vector2 = createVector();
		assertEquals(vector.getLength(), n);
		assertEquals(vector.getOrder(), 1);
		assertEquals(vector.getShape(), new IntCoordinates(n));
		ImmutableMap<IntCoordinates, Double> entries = vector.getCoordinateEntryList();
		for(int i = 0; i < n; i++)
			assertEquals(entries.get(new IntCoordinates(i)).doubleValue(), i);
		assertEquals(vector.size(), n);
		assertEquals(vector, vector2);
	}
	@Test
	public void testAt()
	{
	}
	@Test
	public void testSet()
	{
	}
	@Test
	public void testSetAll()
	{
	}
	@Test
	public void testAdd()
	{
	}
	@Test
	public void testAddAll()
	{
	}
	@Test
	public void testAddInPlace()
	{
	}
	@Test
	public void testMulInPlace()
	{
	}
	@Test
	public void testAddVector()
	{
	}
	@Test
	public void testMul()
	{
	}
	@Test
	public void testOuter()
	{
	}
	@Test
	public void testInner()
	{
	}
	@Test
	public void testMatrixVectorProduct()
	{
	}
}
