package tensorproduct;

import basic.DoubleCompare;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.junit.Test;
import tensorproduct.geometry.CartesianGrid;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class TPShapeFunctionTest
{
	
	@Test
	public void deg1Value()
	{
		final CartesianGrid c = new CartesianGrid(CoordinateVector.fromValues(0, 0), CoordinateVector.fromValues(2,
		                                                                                                         4),
		                                          new IntCoordinates(1, 1));
		TPShapeFunction sf = new TPShapeFunction(c.cells.get(0), 1, 0);
		assertTrue(DoubleCompare.almostEqual(1.0, sf.value(CoordinateVector.fromValues(0, 0))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(2, 0))));
		assertTrue(DoubleCompare.almostEqual(0.5, sf.value(CoordinateVector.fromValues(1, 0))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(0, 4))));
		assertTrue(DoubleCompare.almostEqual(0.25, sf.value(CoordinateVector.fromValues(1, 2))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(2, 4))));
		assertTrue(DoubleCompare.almostEqual(0.5, sf.value(CoordinateVector.fromValues(0, 2))));
		sf = new TPShapeFunction(c.cells.get(0), 1, 1);
		assertTrue(DoubleCompare.almostEqual(1.0, sf.value(CoordinateVector.fromValues(2, 0))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(0, 0))));
		assertTrue(DoubleCompare.almostEqual(0.5, sf.value(CoordinateVector.fromValues(1, 0))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(0, 4))));
		assertTrue(DoubleCompare.almostEqual(0.25, sf.value(CoordinateVector.fromValues(1, 2))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(2, 4))));
		assertTrue(DoubleCompare.almostEqual(0.5, sf.value(CoordinateVector.fromValues(2, 2))));
		sf = new TPShapeFunction(c.cells.get(0), 1, 2);
		assertTrue(DoubleCompare.almostEqual(1.0, sf.value(CoordinateVector.fromValues(0, 4))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(2, 0))));
		assertTrue(DoubleCompare.almostEqual(0.5, sf.value(CoordinateVector.fromValues(1, 4))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(0, 0))));
		assertTrue(DoubleCompare.almostEqual(0.25, sf.value(CoordinateVector.fromValues(1, 2))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(2, 4))));
		assertTrue(DoubleCompare.almostEqual(0.5, sf.value(CoordinateVector.fromValues(0, 2))));
		sf = new TPShapeFunction(c.cells.get(0), 1, 3);
		assertTrue(DoubleCompare.almostEqual(1.0, sf.value(CoordinateVector.fromValues(2, 4))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(2, 0))));
		assertTrue(DoubleCompare.almostEqual(0.5, sf.value(CoordinateVector.fromValues(1, 4))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(0, 4))));
		assertTrue(DoubleCompare.almostEqual(0.25, sf.value(CoordinateVector.fromValues(1, 2))));
		assertTrue(DoubleCompare.almostEqual(0.0, sf.value(CoordinateVector.fromValues(0, 0))));
		assertTrue(DoubleCompare.almostEqual(0.5, sf.value(CoordinateVector.fromValues(2, 2))));
	}
	
	@Test
	public void deg1Grad()
	{
		final CartesianGrid c = new CartesianGrid(CoordinateVector.fromValues(0, 0), CoordinateVector.fromValues(2,
		                                                                                                         4),
		                                          new IntCoordinates(1, 1));
		final TPShapeFunction sf = new TPShapeFunction(c.cells.get(0), 1, 0);
		assertEquals(CoordinateVector.fromValues(0, 0), sf.gradient(CoordinateVector.fromValues(2, 4)));
		assertEquals(CoordinateVector.fromValues(-0.5, 0), sf.gradient(CoordinateVector.fromValues(2, 0)));
		assertEquals(CoordinateVector.fromValues(-0.5, -0.25), sf.gradient(CoordinateVector.fromValues(0, 0)));
		assertEquals(CoordinateVector.fromValues(0, -0.25), sf.gradient(CoordinateVector.fromValues(0, 4)));
	}
	
	@Test
	public void deg1Grad2()
	{
		final CartesianGrid c = new CartesianGrid(CoordinateVector.fromValues(0, 0), CoordinateVector.fromValues(1,
		                                                                                                         1),
		                                          new IntCoordinates(1, 1));
		final TPShapeFunction sf = new TPShapeFunction(c.cells.get(0), 1, 0);
		assertEquals(CoordinateVector.fromValues(-1, -1), sf.gradient(CoordinateVector.fromValues(0, 0)));
		assertEquals(CoordinateVector.fromValues(0, -1), sf.gradient(CoordinateVector.fromValues(0, 1)));
		assertEquals(CoordinateVector.fromValues(0, 0), sf.gradient(CoordinateVector.fromValues(1, 1)));
		assertEquals(CoordinateVector.fromValues(-1, 0), sf.gradient(CoordinateVector.fromValues(1, 0)));
	}
}
