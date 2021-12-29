package mixed;

import basic.DoubleCompare;
import linalg.CoordinateVector;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class MixedValueTest
{
	
	@Test
	public void testMixed()
	{
		final double a = Math.random() * 100;
		final double b = Math.random() * 100;
		final double c = Math.random() * 100;
		final double d = Math.random() - 1;
		final double e = Math.random() - 1;
		final double f = Math.random() - 1;
		final PressureValue p21 = new PressureValue(d, 2);
		final PressureValue p22 = new PressureValue(a, 2);
		final PressureValue p31 = new PressureValue(d, 3);
		final PressureValue p32 = new PressureValue(a, 3);
		final VelocityValue v21 = new VelocityValue(CoordinateVector.fromValues(d, e));
		final VelocityValue v22 = new VelocityValue(CoordinateVector.fromValues(a, b));
		final VelocityValue v31 = new VelocityValue(CoordinateVector.fromValues(d, e, f));
		final VelocityValue v32 = new VelocityValue(CoordinateVector.fromValues(a, b, c));
		MixedValue m1 = p21.add(v21);
		MixedValue m2 = p31.add(v31);
		final MixedValue m3 = p32.add(v31);
		m1 = m1.mul(2.0);
		m1 = m1.add(v22);
		m2 = m2.add(v32);
		m2.add(c, 2);
		m1.addPressure(p22.getPressure());
		m2 = m2.add(m3);
		assertTrue(DoubleCompare.almostEqual(m1.at(0), d * 2 + a));
		assertTrue(DoubleCompare.almostEqual(m1.at(1), d * 2 + a));
		assertTrue(DoubleCompare.almostEqual(m1.at(2), e * 2 + b));
		assertTrue(DoubleCompare.almostEqual(m2.at(0), d + a));
		assertTrue(Math.abs(m2.at(1) - 2 * d - a) < 1e-12);
		assertTrue(Math.abs(m2.at(2) - 2 * e - b - c) < 1e-12);
		assertTrue(Math.abs(m2.at(3) - 2 * f - c) < 1e-12);
		assertTrue(DoubleCompare.almostEqual(m3.at(0), a));
		assertTrue(DoubleCompare.almostEqual(m3.at(1), d));
		assertTrue(DoubleCompare.almostEqual(m3.at(2), e));
		assertTrue(DoubleCompare.almostEqual(m3.at(3), f));
	}
	
	@Test
	public void testPure()
	{
		final double a = Math.random() * 100;
		final double b = Math.random() * 100;
		final double c = Math.random() * 100;
		final double d = Math.random() - 1;
		final double e = Math.random() - 1;
		final double f = Math.random() - 1;
		final PressureValue p21 = new PressureValue(d, 2);
		final PressureValue p22 = new PressureValue(0, 2);
		final PressureValue p31 = new PressureValue(d, 3);
		final PressureValue p32 = new PressureValue(0, 3);
		final VelocityValue v21 = new VelocityValue(CoordinateVector.fromValues(d, e));
		final VelocityValue v22 = new VelocityValue(CoordinateVector.fromValues(a, 0));
		VelocityValue v31 = new VelocityValue(CoordinateVector.fromValues(d, e, f));
		final VelocityValue v32 = new VelocityValue(CoordinateVector.fromValues(a, b, c));
		p21.add(d, 0);
		p31.add(d, 0);
		p21.addPressure(d);
		p22.setPressure(d);
		p22.set(a, 0);
		p31.addPressure(d);
		p32.setPressure(d);
		p32.set(a, 0);
		v21.add(d, 1);
		v22.setVelocity(CoordinateVector.fromValues(a, b));
		final MixedValue v33 = v32.add(v31);
		v31 = v31.mul(3);
		assertEquals(p21.getLength(), 3);
		assertEquals(p22.getLength(), 3);
		assertEquals(p31.getLength(), 4);
		assertEquals(p32.getLength(), 4);
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
		assertTrue(DoubleCompare.almostEqual(p21.at(0), 3 * d));
		assertTrue(DoubleCompare.almostEqual(p22.at(0), a));
		assertTrue(DoubleCompare.almostEqual(p31.at(0), 3 * d));
		assertTrue(DoubleCompare.almostEqual(p32.at(0), a));
		assertTrue(DoubleCompare.almostEqual(v21.at(1), 2 * d));
		assertTrue(DoubleCompare.almostEqual(v21.at(2), e));
		assertTrue(DoubleCompare.almostEqual(v22.at(1), a));
		assertTrue(DoubleCompare.almostEqual(v22.at(2), b));
		assertTrue(DoubleCompare.almostEqual(v31.at(1), d * 3));
		assertTrue(DoubleCompare.almostEqual(v31.at(2), e * 3));
		assertTrue(DoubleCompare.almostEqual(v31.at(3), f * 3));
		assertTrue(DoubleCompare.almostEqual(v32.at(1), a));
		assertTrue(DoubleCompare.almostEqual(v32.at(2), b));
		assertTrue(DoubleCompare.almostEqual(v32.at(3), c));
		assertTrue(DoubleCompare.almostEqual(v33.at(1), a + d));
		assertTrue(DoubleCompare.almostEqual(v33.at(2), b + e));
		assertTrue(DoubleCompare.almostEqual(v33.at(3), c + f));
	}
}
