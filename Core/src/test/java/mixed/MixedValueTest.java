package mixed;

import linalg.CoordinateVector;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class MixedValueTest
{
	
	@Test
	public void testMixed()
	{
		double a = Math.random() * 100;
		double b = Math.random() * 100;
		double c = Math.random() * 100;
		double d = Math.random() - 1;
		double e = Math.random() - 1;
		double f = Math.random() - 1;
		PressureValue p1 = new PressureValue(d);
		PressureValue p2 = new PressureValue(a);
		VelocityValue v21 = new VelocityValue(CoordinateVector.fromValues(d, e));
		VelocityValue v22 = new VelocityValue(CoordinateVector.fromValues(a,b));
		VelocityValue v31 = new VelocityValue(CoordinateVector.fromValues(d,e,f));
		VelocityValue v32 = new VelocityValue(CoordinateVector.fromValues(a,b,c));
		MixedValue m1 = p1.add(v21);
		MixedValue m2 = p1.add(v31);
		MixedValue m3 = p2.add(v31);
		m1 = m1.mul(2.0);
		m1 = m1.add(v22);
		m2 = m2.add(v32);
		m2.add(c,2);
		m1.addPressure(p2.getPressure());
		m2 = m2.add(m3);
		assertEquals(m1.at(0),d*2+a);
		assertEquals(m1.at(1),d*2+a);
		assertEquals(m1.at(2),e*2+b);
		assertEquals(m2.at(0),d+a);
		assertTrue(Math.abs(m2.at(1)-2*d-a)<1e-12);
		assertTrue(Math.abs(m2.at(2)-2*e-b-c)<1e-12);
		assertTrue(Math.abs(m2.at(3)-2*f-c)<1e-12);
		assertEquals(m3.at(0),a);
		assertEquals(m3.at(1),d);
		assertEquals(m3.at(2),e);
		assertEquals(m3.at(3),f);
		
	}
	@Test
	public void testPure() {
		double a = Math.random()*100;
		double b = Math.random()*100;
		double c = Math.random()*100;
		double d = Math.random()-1;
		double e = Math.random()-1;
		double f = Math.random()-1;
		PressureValue p1 = new PressureValue(d);
		PressureValue p2 = new PressureValue(0);
		VelocityValue v21 = new VelocityValue(CoordinateVector.fromValues(d, e));
		VelocityValue v22 = new VelocityValue(CoordinateVector.fromValues(a,0));
		VelocityValue v31 = new VelocityValue(CoordinateVector.fromValues(d,e,f));
		VelocityValue v32 = new VelocityValue(CoordinateVector.fromValues(a,b,c));
		p1.add(d,0);
		p1.addPressure(d);
		p2.setPressure(d);
		p2.set(a,0);
		v21.add(d,1);
		v22.setVelocity(CoordinateVector.fromValues(a,b));
		MixedValue v33 = v32.add(v31);
		v31 = v31.mul(3);
		assertEquals(p1.getLength(),4);
		assertEquals(p2.getLength(), 4);
		assertEquals(v21.getLength(), 3);
		assertEquals(v22.getLength(), 3);
		assertEquals(v31.getLength(), 4);
		assertEquals(v32.getLength(), 4);
		assertEquals(v32.getLength(), 4);
		assertEquals(v21.getDomainDimension(), 2);
		assertEquals(v22.getDomainDimension(), 2);
		assertEquals(v31.getDomainDimension(), 3);
		assertEquals(v32.getDomainDimension(), 3);
		assertEquals(v33.getDomainDimension(), 3);
		assertEquals(p1.at(0),3*d);
		assertEquals(p2.at(0),a);
		assertEquals(v21.at(1),2*d);
		assertEquals(v21.at(2),e);
		assertEquals(v22.at(1),a);
		assertEquals(v22.at(2),b);
		assertEquals(v31.at(1),d*3);
		assertEquals(v31.at(2),e*3);
		assertEquals(v31.at(3),f*3);
		assertEquals(v32.at(1),a);
		assertEquals(v32.at(2),b);
		assertEquals(v32.at(3),c);
		assertEquals(v33.at(1),a+d);
		assertEquals(v33.at(2),b+e);
		assertEquals(v33.at(3),c+f);
	}
	
	
}