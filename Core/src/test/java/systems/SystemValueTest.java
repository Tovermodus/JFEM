package systems;

import linalg.CoordinateVector;
import mixed.MixedValue;
import mixed.PressureValue;
import mixed.VelocityValue;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class SystemValueTest
{
	
	@Test
	public void testGetComponentZeros()
	{
		int [] ends = new int[]{1,4,6,7};
		SystemValue v = new SystemValue(ends);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), new CoordinateVector(1));
		assertEquals(v.getComponent(1), new CoordinateVector(3));
		assertEquals(v.getComponent(2), new CoordinateVector(2));
		assertEquals(v.getComponent(3), new CoordinateVector(1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(5));
	}
	@Test
	public void testGetComponent()
	{
		int [] ends = new int[]{1,4,6,7};
		double [] values = new double[]{0.1,0.5,0.4,0.3,0.8,0.7,1.0};
		SystemValue v = new SystemValue(ends, values);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), CoordinateVector.fromValues(0.1));
		assertEquals(v.getComponent(1), CoordinateVector.fromValues(0.5,0.4,0.3));
		assertEquals(v.getComponent(2), CoordinateVector.fromValues(0.8,0.7));
		assertEquals(v.getComponent(3), CoordinateVector.fromValues(1.0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(5));
		
	}
	@Test
	public void testGetComponentParameters()
	{
		int [] ends = new int[]{1,4,6,7};
		SystemParameters.createInstance(ends);
		double [] values = new double[]{0.1,0.5,0.4,0.3,0.8,0.7,1.0};
		SystemValue v = new SystemValue( values);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), CoordinateVector.fromValues(0.1));
		assertEquals(v.getComponent(1), CoordinateVector.fromValues(0.5,0.4,0.3));
		assertEquals(v.getComponent(2), CoordinateVector.fromValues(0.8,0.7));
		assertEquals(v.getComponent(3), CoordinateVector.fromValues(1.0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(5));
		SystemParameters.deleteInstance();
		
	}
	@Test
	public void testSetComponentZeros() {
		int [] ends = new int[]{1,4,6,7};
		SystemValue v = new SystemValue(ends);
		CoordinateVector c1 = CoordinateVector.fromValues(0.999,0.99);
		CoordinateVector c2 = CoordinateVector.fromValues(777,7777,77777);
		v.setComponent(c1,2);
		v.setComponent(c2,1);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), new CoordinateVector(1));
		assertEquals(v.getComponent(1), c2);
		assertEquals(v.getComponent(2), c1);
		assertEquals(v.getComponent(3), new CoordinateVector(1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(5));
	
	}
	@Test
	public void testSetComponent() {
		int [] ends = new int[]{1,4,6,7};
		SystemParameters.createInstance(ends);
		double [] values = new double[]{0.1,0.5,0.4,0.3,0.8,0.7,1.0};
		SystemValue v = new SystemValue( values);
		CoordinateVector c1 = CoordinateVector.fromValues(0.999,0.99);
		CoordinateVector c2 = CoordinateVector.fromValues(777,7777,77777);
		v.setComponent(c1,2);
		v.setComponent(c2,1);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), CoordinateVector.fromValues(0.1));
		assertEquals(v.getComponent(1), c2);
		assertEquals(v.getComponent(2), c1);
		assertEquals(v.getComponent(3), CoordinateVector.fromValues(1.0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(5));
		SystemParameters.deleteInstance();
	
	}
	@Test
	public void testSetGetComponent() {
		int [] ends = new int[]{1,4,6,7};
		double [] values = new double[]{0.1,0.5,0.4,0.3,0.8,0.7,1.0};
		SystemValue v = new SystemValue(ends, values);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), CoordinateVector.fromValues(0.1));
		assertEquals(v.getComponent(1), CoordinateVector.fromValues(0.5,0.4,0.3));
		assertEquals(v.getComponent(2), CoordinateVector.fromValues(0.8,0.7));
		assertEquals(v.getComponent(3), CoordinateVector.fromValues(1.0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(5));
		CoordinateVector c1 = CoordinateVector.fromValues(0.999,0.99);
		CoordinateVector c2 = CoordinateVector.fromValues(777,7777,77777);
		v.setComponent(c1,2);
		v.setComponent(c2,1);
		assertEquals(v.getNumberOfComponents(), 4);
		assertEquals(v.getComponent(0), CoordinateVector.fromValues(0.1));
		assertEquals(v.getComponent(1), c2);
		assertEquals(v.getComponent(2), c1);
		assertEquals(v.getComponent(3), CoordinateVector.fromValues(1.0));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(-1));
		assertThrows(IllegalArgumentException.class, ()->v.getComponent(5));
	
	}
	
	
}