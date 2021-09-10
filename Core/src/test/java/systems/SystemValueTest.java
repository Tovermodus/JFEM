package systems;

import linalg.CoordinateVector;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

public class SystemValueTest
{
	
	@Test
	public void testGetComponentZeros()
	{
		final int[] ends = new int[]{1, 4, 6, 7};
		final SystemValue v = new SystemValue(ends);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), new CoordinateVector(1));
		assertEquals(v.getComponent(1), new CoordinateVector(3));
		assertEquals(v.getComponent(2), new CoordinateVector(2));
		assertEquals(v.getComponent(3), new CoordinateVector(1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(4));
	}
	
	@Test
	public void testGetComponent()
	{
		final int[] ends = new int[]{1, 4, 6, 7};
		final double[] values = new double[]{0.1, 0.5, 0.4, 0.3, 0.8, 0.7, 1.0};
		final SystemValue v = new SystemValue(ends, values);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), CoordinateVector.fromValues(0.1));
		assertEquals(v.getComponent(1), CoordinateVector.fromValues(0.5, 0.4, 0.3));
		assertEquals(v.getComponent(2), CoordinateVector.fromValues(0.8, 0.7));
		assertEquals(v.getComponent(3), CoordinateVector.fromValues(1.0));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(5));
	}
	
	@Test
	public void testGetComponentParameters()
	{
		final int[] ends = new int[]{1, 4, 6, 7};
		SystemParameters.createInstance(ends);
		final double[] values = new double[]{0.1, 0.5, 0.4, 0.3, 0.8, 0.7, 1.0};
		final SystemValue v = new SystemValue(values);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), CoordinateVector.fromValues(0.1));
		assertEquals(v.getComponent(1), CoordinateVector.fromValues(0.5, 0.4, 0.3));
		assertEquals(v.getComponent(2), CoordinateVector.fromValues(0.8, 0.7));
		assertEquals(v.getComponent(3), CoordinateVector.fromValues(1.0));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(5));
		SystemParameters.deleteInstance();
	}
	
	@Test
	public void testSetComponentZeros()
	{
		final int[] ends = new int[]{1, 4, 6, 7};
		final SystemValue v = new SystemValue(ends);
		final CoordinateVector c1 = CoordinateVector.fromValues(0.999, 0.99);
		final CoordinateVector c2 = CoordinateVector.fromValues(777, 7777, 77777);
		v.setComponent(c1, 2);
		v.setComponent(c2, 1);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), new CoordinateVector(1));
		assertEquals(v.getComponent(1), c2);
		assertEquals(v.getComponent(2), c1);
		assertEquals(v.getComponent(3), new CoordinateVector(1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(5));
	}
	
	@Test
	public void testSetComponent()
	{
		final int[] ends = new int[]{1, 4, 6, 7};
		SystemParameters.createInstance(ends);
		final double[] values = new double[]{0.1, 0.5, 0.4, 0.3, 0.8, 0.7, 1.0};
		final SystemValue v = new SystemValue(values);
		final CoordinateVector c1 = CoordinateVector.fromValues(0.999, 0.99);
		final CoordinateVector c2 = CoordinateVector.fromValues(777, 7777, 77777);
		v.setComponent(c1, 2);
		v.setComponent(c2, 1);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), CoordinateVector.fromValues(0.1));
		assertEquals(v.getComponent(1), c2);
		assertEquals(v.getComponent(2), c1);
		assertEquals(v.getComponent(3), CoordinateVector.fromValues(1.0));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(5));
		SystemParameters.deleteInstance();
	}
	
	@Test
	public void testSetGetComponent()
	{
		final int[] ends = new int[]{1, 4, 6, 7};
		final double[] values = new double[]{0.1, 0.5, 0.4, 0.3, 0.8, 0.7, 1.0};
		final SystemValue v = new SystemValue(ends, values);
		v.setComponent(-0.1, 0);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), CoordinateVector.fromValues(-0.1));
		assertEquals(v.getComponent(1), CoordinateVector.fromValues(0.5, 0.4, 0.3));
		assertEquals(v.getComponent(2), CoordinateVector.fromValues(0.8, 0.7));
		assertEquals(v.getComponent(3), CoordinateVector.fromValues(1.0));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(5));
		final CoordinateVector c1 = CoordinateVector.fromValues(0.999, 0.99);
		final CoordinateVector c2 = CoordinateVector.fromValues(777, 7777, 77777);
		v.setComponent(c1, 2);
		v.setComponent(c2, 1);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), CoordinateVector.fromValues(-0.1));
		assertEquals(v.getComponent(1), c2);
		assertEquals(v.getComponent(2), c1);
		assertEquals(v.getComponent(3), CoordinateVector.fromValues(1.0));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(5));
	}
}
