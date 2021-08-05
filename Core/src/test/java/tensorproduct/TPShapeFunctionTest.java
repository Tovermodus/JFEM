package tensorproduct;

import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.junit.jupiter.api.Test;
import tensorproduct.geometry.CartesianGrid;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class TPShapeFunctionTest
{
	
	@Test
	public void deg1Value()
	{
		CartesianGrid c = new CartesianGrid(CoordinateVector.fromValues(0,0), CoordinateVector.fromValues(2,
			4), new IntCoordinates(1,1));
		TPShapeFunction sf = new TPShapeFunction(c.cells.get(0), 1, 0);
		assertEquals(1.0, sf.value(CoordinateVector.fromValues(0,0)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(2,0)).doubleValue());
		assertEquals(0.5, sf.value(CoordinateVector.fromValues(1,0)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(0,4)).doubleValue());
		assertEquals(0.25, sf.value(CoordinateVector.fromValues(1,2)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(2,4)).doubleValue());
		assertEquals(0.5, sf.value(CoordinateVector.fromValues(0,2)).doubleValue());
		sf = new TPShapeFunction(c.cells.get(0), 1, 1);
		assertEquals(1.0, sf.value(CoordinateVector.fromValues(2,0)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(0,0)).doubleValue());
		assertEquals(0.5, sf.value(CoordinateVector.fromValues(1,0)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(0,4)).doubleValue());
		assertEquals(0.25, sf.value(CoordinateVector.fromValues(1,2)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(2,4)).doubleValue());
		assertEquals(0.5, sf.value(CoordinateVector.fromValues(2,2)).doubleValue());
		sf = new TPShapeFunction(c.cells.get(0), 1, 2);
		assertEquals(1.0, sf.value(CoordinateVector.fromValues(0,4)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(2,0)).doubleValue());
		assertEquals(0.5, sf.value(CoordinateVector.fromValues(1,4)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(0,0)).doubleValue());
		assertEquals(0.25, sf.value(CoordinateVector.fromValues(1,2)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(2,4)).doubleValue());
		assertEquals(0.5, sf.value(CoordinateVector.fromValues(0,2)).doubleValue());
		sf = new TPShapeFunction(c.cells.get(0), 1, 3);
		assertEquals(1.0, sf.value(CoordinateVector.fromValues(2,4)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(2,0)).doubleValue());
		assertEquals(0.5, sf.value(CoordinateVector.fromValues(1,4)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(0,4)).doubleValue());
		assertEquals(0.25, sf.value(CoordinateVector.fromValues(1,2)).doubleValue());
		assertEquals(0.0, sf.value(CoordinateVector.fromValues(0,0)).doubleValue());
		assertEquals(0.5, sf.value(CoordinateVector.fromValues(2,2)).doubleValue());
		
	}
	
	@Test
	public void deg1Grad()
	{
		CartesianGrid c = new CartesianGrid(CoordinateVector.fromValues(0,0), CoordinateVector.fromValues(2,
			4), new IntCoordinates(1,1));
		TPShapeFunction sf = new TPShapeFunction(c.cells.get(0), 1, 0);
		assertEquals(CoordinateVector.fromValues(0,0), sf.gradient(CoordinateVector.fromValues(2,4)));
		assertEquals(CoordinateVector.fromValues(-0.5,0), sf.gradient(CoordinateVector.fromValues(2,0)));
		assertEquals(CoordinateVector.fromValues(-0.5,-0.25), sf.gradient(CoordinateVector.fromValues(0,0)));
		assertEquals(CoordinateVector.fromValues(0,-0.25), sf.gradient(CoordinateVector.fromValues(0,4)));
		
	}
	@Test
	public void deg1Grad2()
	{
		CartesianGrid c = new CartesianGrid(CoordinateVector.fromValues(0,0), CoordinateVector.fromValues(1,
			1), new IntCoordinates(1,1));
		TPShapeFunction sf = new TPShapeFunction(c.cells.get(0), 1, 0);
		assertEquals(CoordinateVector.fromValues(-1,-1), sf.gradient(CoordinateVector.fromValues(0,0)));
		assertEquals(CoordinateVector.fromValues(0,-1), sf.gradient(CoordinateVector.fromValues(0,1)));
		assertEquals(CoordinateVector.fromValues(0,0), sf.gradient(CoordinateVector.fromValues(1,1)));
		assertEquals(CoordinateVector.fromValues(-1,0), sf.gradient(CoordinateVector.fromValues(1,0)));
		
	}
}
