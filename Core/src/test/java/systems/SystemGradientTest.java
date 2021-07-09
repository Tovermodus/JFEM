package systems;

import linalg.CoordinateMatrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class SystemGradientTest
{
	public SystemGradient createSystemGradient()
	{
		int [] ends = new int[]{1,3};
		SystemGradient v = new SystemGradient(ends);
		v.set(0.1,0,0);
		v.set(0.3,0,1);
		v.set(0.4,0,2);
		v.set(0.6,1,0);
		v.set(0.9,1,1);
		v.set(1.0,1,2);
		v.set(0.7,2,0);
		v.set(1.1,2,1);
		v.set(1.2,2,2);
		return v;
	}
	@Test
	public void testGetComponentZeros()
	{
		int [] ends = new int[]{1,2,4};
		SystemGradient v = new SystemGradient(ends);
		assertEquals(v.getNumberOfComponents(), 3);
		assertEquals(v.getComponent(0,0), new CoordinateMatrix(1,1));
		assertEquals(v.getComponent(0,1), new CoordinateMatrix(1,1));
		assertEquals(v.getComponent(0,2), new CoordinateMatrix(1,2));
		assertEquals(v.getComponent(1,0), new CoordinateMatrix(1,1));
		assertEquals(v.getComponent(1,1), new CoordinateMatrix(1,1));
		assertEquals(v.getComponent(1,2), new CoordinateMatrix(1,2));
		assertEquals(v.getComponent(2,0), new CoordinateMatrix(2,1));
		assertEquals(v.getComponent(2,1), new CoordinateMatrix(2,1));
		assertEquals(v.getComponent(2,2), new CoordinateMatrix(2,2));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1,2));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(2,-1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(3,0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(1,3));
	}
	@Test
	public void testGetComponent()
	{
		SystemGradient v = createSystemGradient();
		assertEquals(v.getNumberOfComponents(), 2);
		assertEquals(v.getComponent(0,0), CoordinateMatrix.fromValues(1,1,0.1));
		assertEquals(v.getComponent(0,1), CoordinateMatrix.fromValues(1,2,0.3,0.4));
		assertEquals(v.getComponent(1,0), CoordinateMatrix.fromValues(2,1,0.6,0.7));
		assertEquals(v.getComponent(1,1), CoordinateMatrix.fromValues(2,2,0.9,1.0,1.1,1.2));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1,0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(5,1));
		
	}
	@Test
	public void testGetComponentZerosParameters()
	{
		int [] ends = new int[]{1,2,4};
		SystemParameters.createInstance(ends);
		SystemGradient v = new SystemGradient();
		assertEquals(v.getNumberOfComponents(), 3);
		assertEquals(v.getComponent(0,0), new CoordinateMatrix(1,1));
		assertEquals(v.getComponent(0,1), new CoordinateMatrix(1,1));
		assertEquals(v.getComponent(0,2), new CoordinateMatrix(1,2));
		assertEquals(v.getComponent(1,0), new CoordinateMatrix(1,1));
		assertEquals(v.getComponent(1,1), new CoordinateMatrix(1,1));
		assertEquals(v.getComponent(1,2), new CoordinateMatrix(1,2));
		assertEquals(v.getComponent(2,0), new CoordinateMatrix(2,1));
		assertEquals(v.getComponent(2,1), new CoordinateMatrix(2,1));
		assertEquals(v.getComponent(2,2), new CoordinateMatrix(2,2));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1,2));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(2,-1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(3,0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(1,3));
		SystemParameters.deleteInstance();
		
	}
	@Test
	public void testSetComponentZeros() {
		int [] ends = new int[]{1,2,4};
		SystemParameters.createInstance(ends);
		SystemGradient v = new SystemGradient();
		CoordinateMatrix c2 = CoordinateMatrix.fromValues(2,2,9.0,0.0,0.0,9.0);
		CoordinateMatrix c1 = CoordinateMatrix.fromValues(1,1,-1000.0);
		v.setComponent(c1, 1,1);
		v.setComponent(c2,2,2);
		assertEquals(v.getNumberOfComponents(), 3);
		assertEquals(v.getComponent(0,0), new CoordinateMatrix(1,1));
		assertEquals(v.getComponent(0,1), new CoordinateMatrix(1,1));
		assertEquals(v.getComponent(0,2), new CoordinateMatrix(1,2));
		assertEquals(v.getComponent(1,0), new CoordinateMatrix(1,1));
		assertTrue(v.getComponent(1,1).almostEqual(c1));
		assertEquals(v.getComponent(1,1), c1);
		assertEquals(v.getComponent(1,2), new CoordinateMatrix(1,2));
		assertEquals(v.getComponent(2,0), new CoordinateMatrix(2,1));
		assertEquals(v.getComponent(2,1), new CoordinateMatrix(2,1));
		assertEquals(v.getComponent(2,2), c2);
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1,2));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(2,-1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(3,0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(1,3));
		SystemParameters.deleteInstance();
	
	}
	@Test
	public void testSetComponent() {
		SystemGradient v = createSystemGradient();
		CoordinateMatrix c1 = CoordinateMatrix.fromValues(2,2,9.0,0.0,0.0,9.0);
		v.setComponent(c1,1,1);
		assertEquals(v.getNumberOfComponents(), 2);
		assertEquals(v.getComponent(0,0), CoordinateMatrix.fromValues(1,1,0.1));
		assertEquals(v.getComponent(0,1), CoordinateMatrix.fromValues(1,2,0.3,0.4));
		assertEquals(v.getComponent(1,0), CoordinateMatrix.fromValues(2,1,0.6,0.7));
		assertEquals(v.getComponent(1,1), c1);
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1,0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(5,1));
	
	}
	@Test
	public void testSetGetComponent() {
		SystemGradient v = createSystemGradient();
		CoordinateMatrix c1 = CoordinateMatrix.fromValues(1,2,9.0,0.0);
		assertEquals(v.getNumberOfComponents(), 2);
		assertEquals(v.getComponent(0,0), CoordinateMatrix.fromValues(1,1,0.1));
		assertEquals(v.getComponent(0,1), CoordinateMatrix.fromValues(1,2,0.3,0.4));
		assertEquals(v.getComponent(1,0), CoordinateMatrix.fromValues(2,1,0.6,0.7));
		assertEquals(v.getComponent(1,1), CoordinateMatrix.fromValues(2,2,0.9,1.0,1.1,1.2));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1,0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(5,1));
		v.setComponent(c1,0,1);
		assertEquals(v.getComponent(0,0), CoordinateMatrix.fromValues(1,1,0.1));
		assertEquals(v.getComponent(0,1), c1);
		assertEquals(v.getComponent(1,0), CoordinateMatrix.fromValues(2,1,0.6,0.7));
		assertEquals(v.getComponent(1,1), CoordinateMatrix.fromValues(2,2,0.9,1.0,1.1,1.2));
	
	}
	
	
}