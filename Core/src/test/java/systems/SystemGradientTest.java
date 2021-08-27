package systems;

import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

public class SystemGradientTest
{
	public static SystemGradient createSystemGradient()
	{
		final int[] ends = new int[]{1, 3};
		final SystemGradient v = new SystemGradient(ends, 3);
		v.set(0.1, 0, 0);
		v.set(0.3, 0, 1);
		v.set(0.4, 0, 2);
		v.set(0.6, 1, 0);
		v.set(0.9, 1, 1);
		v.set(1.0, 1, 2);
		v.set(0.7, 2, 0);
		v.set(1.1, 2, 1);
		v.set(1.2, 2, 2);
		return v;
	}
	
	@Test
	public void testGetComponentZeros()
	{
		final int[] ends = new int[]{1, 2, 4};
		final SystemGradient v = new SystemGradient(ends, 3);
		assertEquals(v.getNumberOfComponents(), 3);
		assertEquals(v.getComponent(0), new CoordinateDenseMatrix(1, 3));
		assertEquals(v.getComponent(1), new CoordinateDenseMatrix(1, 3));
		assertEquals(v.getComponent(2), new CoordinateDenseMatrix(2, 3));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(3));
	}
	
	@Test
	public void testGetComponent()
	{
		final SystemGradient v = createSystemGradient();
		assertEquals(v.getNumberOfComponents(), 2);
		assertEquals(v.getComponent(0), CoordinateDenseMatrix.fromValues(1, 3, 0.1, 0.3, 0.4));
		assertEquals(v.getComponent(1), CoordinateDenseMatrix.fromValues(2, 3, 0.6, 0.9, 1.0, 0.7, 1.1, 1.2));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(3));
	}
	
	@Test
	public void testGetComponentZerosParameters()
	{
		final int[] ends = new int[]{1, 2, 4};
		SystemParameters.createInstance(ends);
		final SystemGradient v = new SystemGradient(3);
		assertEquals(v.getNumberOfComponents(), 3);
		assertEquals(v.getComponent(0), new CoordinateDenseMatrix(1, 3));
		assertEquals(v.getComponent(1), new CoordinateDenseMatrix(1, 3));
		assertEquals(v.getComponent(2), new CoordinateDenseMatrix(2, 3));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(3));
		SystemParameters.deleteInstance();
	}
	
	@Test
	public void testSetComponentZeros()
	{
		final int[] ends = new int[]{1, 2, 4};
		SystemParameters.createInstance(ends);
		final SystemGradient v = new SystemGradient(3);
		final CoordinateMatrix c2 = CoordinateDenseMatrix.fromValues(2, 3, 9.0, 0.0, 0.0, 9.0, 9.0, 0.0);
		final CoordinateVector c1 = CoordinateVector.fromValues(-1000.0, -2000.0, -3000.0);
		v.setComponent(c1, 1);
		v.setComponent(c2, 2);
		assertEquals(v.getNumberOfComponents(), 3);
		assertEquals(v.getComponent(0), new CoordinateDenseMatrix(1, 3));
		assertEquals(v.getComponent(1), c1.asCoordinateMatrix().transpose());
		assertEquals(v.getComponent(2), c2);
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(3));
		SystemParameters.deleteInstance();
	}
	
	@Test
	public void testSetComponent()
	{
		final SystemGradient v = createSystemGradient();
		final CoordinateMatrix c2 = CoordinateDenseMatrix.fromValues(2, 3, 9.0, 0.0, 0.0, 9.0, 9.0, 0.0);
		final CoordinateMatrix c1 = CoordinateVector
			.fromValues(-1000.0, -2000.0, -3000.0)
			.asCoordinateMatrix()
			.transpose();
		v.setComponent(c1, 0);
		v.setComponent(c2, 1);
		assertEquals(v.getNumberOfComponents(), 2);
		assertEquals(v.getComponent(0), c1);
		assertEquals(v.getComponent(1), c2);
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(3));
		SystemParameters.deleteInstance();
	}
	
	@Test
	public void testSetGetComponent()
	{
		final SystemGradient v = createSystemGradient();
		assertEquals(v.getComponent(0), CoordinateDenseMatrix.fromValues(1, 3, 0.1, 0.3, 0.4));
		assertEquals(v.getComponent(1), CoordinateDenseMatrix.fromValues(2, 3, 0.6, 0.9, 1.0, 0.7, 1.1, 1.2));
		assertEquals(v.getNumberOfComponents(), 2);
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(3));
		final CoordinateMatrix c2 = CoordinateDenseMatrix.fromValues(2, 3, 9.0, 0.0, 0.0, 9.0, 9.0, 0.0);
		final CoordinateVector c1 = CoordinateVector.fromValues(-1000.0, -2000.0, -3000.0);
		assertEquals(v.getNumberOfComponents(), 2);
		v.setComponent(c1, 0);
		v.setComponent(c2, 1);
		assertEquals(v.getComponent(0), c1.asCoordinateMatrix().transpose());
		assertEquals(v.getComponent(1), c2);
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, () -> v.getComponent(3));
	}
}
