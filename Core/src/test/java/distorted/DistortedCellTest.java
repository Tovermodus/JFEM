package distorted;

import distorted.geometry.DistortedCell;
import linalg.CoordinateVector;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

public class DistortedCellTest
{
	@Test
	public void testOrdering2D()
	{
		CoordinateVector [] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0,0);
		vertices[1] = CoordinateVector.fromValues(1,1);
		vertices[2] = CoordinateVector.fromValues(0,1);
		vertices[3] = CoordinateVector.fromValues(1,0);
		assertThrows(IllegalArgumentException.class, ()->new DistortedCell(vertices));
		
		vertices[1] = CoordinateVector.fromValues(0,0);
		vertices[2] = CoordinateVector.fromValues(1,0);
		vertices[3] = CoordinateVector.fromValues(1,1);
		vertices[0] = CoordinateVector.fromValues(0,1);
		assertThrows(IllegalArgumentException.class, ()->new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0,0);
		vertices[1] = CoordinateVector.fromValues(0.5,1.5);
		vertices[2] = CoordinateVector.fromValues(2,2);
		vertices[3] = CoordinateVector.fromValues(0,2);
		assertThrows(IllegalArgumentException.class, ()->new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0,0);
		vertices[1] = CoordinateVector.fromValues(0,2);
		vertices[2] = CoordinateVector.fromValues(2,2);
		vertices[3] = CoordinateVector.fromValues(2,0);
		assertThrows(IllegalArgumentException.class, ()->new DistortedCell(vertices));
	}
	
	@Test
	public void testOrdering3D()
	{
		
		CoordinateVector [] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0,0,0);
		vertices[1] = CoordinateVector.fromValues(-1,0,0);
		vertices[2] = CoordinateVector.fromValues(1,1,0);
		vertices[3] = CoordinateVector.fromValues(0,1,0);
		vertices[4] = CoordinateVector.fromValues(0,0,1);
		vertices[5] = CoordinateVector.fromValues(1,0,1);
		vertices[6] = CoordinateVector.fromValues(1,1,1);
		vertices[7] = CoordinateVector.fromValues(0,1,1);
		assertThrows(IllegalArgumentException.class, ()->new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0,0,0);
		vertices[1] = CoordinateVector.fromValues(1,0,2);
		vertices[2] = CoordinateVector.fromValues(1,1,0);
		vertices[3] = CoordinateVector.fromValues(0,1,0);
		vertices[4] = CoordinateVector.fromValues(0,0,1);
		vertices[5] = CoordinateVector.fromValues(1,0,1);
		vertices[6] = CoordinateVector.fromValues(1,1,1);
		vertices[7] = CoordinateVector.fromValues(0,1,1);
		assertThrows(IllegalArgumentException.class, ()->new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0,0,0);
		vertices[1] = CoordinateVector.fromValues(1,0,0);
		vertices[2] = CoordinateVector.fromValues(1,1,0);
		vertices[3] = CoordinateVector.fromValues(0,1,0);
		vertices[4] = CoordinateVector.fromValues(0,0,1);
		vertices[5] = CoordinateVector.fromValues(0,1,1);
		vertices[6] = CoordinateVector.fromValues(1,1,1);
		vertices[7] = CoordinateVector.fromValues(1,0,1);
		assertThrows(IllegalArgumentException.class, ()->new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0,0,0);
		vertices[1] = CoordinateVector.fromValues(1,0,0);
		vertices[2] = CoordinateVector.fromValues(1,1,0);
		vertices[3] = CoordinateVector.fromValues(0,1,0);
		vertices[4] = CoordinateVector.fromValues(0,0,1);
		vertices[5] = CoordinateVector.fromValues(1,0,1);
		vertices[6] = CoordinateVector.fromValues(1,1,-1);
		vertices[7] = CoordinateVector.fromValues(0,1,1);
		assertThrows(IllegalArgumentException.class, ()->new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0,0,0);
		vertices[1] = CoordinateVector.fromValues(1,0,0);
		vertices[2] = CoordinateVector.fromValues(1,1,0);
		vertices[3] = CoordinateVector.fromValues(0,1,0);
		vertices[4] = CoordinateVector.fromValues(1,0,1);
		vertices[5] = CoordinateVector.fromValues(1,1,1);
		vertices[6] = CoordinateVector.fromValues(0,1,1);
		vertices[7] = CoordinateVector.fromValues(0,0,1);
		assertThrows(IllegalArgumentException.class, ()->new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0,0,0);
		vertices[1] = CoordinateVector.fromValues(1,0,0);
		vertices[2] = CoordinateVector.fromValues(1,1,1);
		vertices[3] = CoordinateVector.fromValues(0,1,1);
		vertices[4] = CoordinateVector.fromValues(0,0,1);
		vertices[5] = CoordinateVector.fromValues(1,0,1);
		vertices[6] = CoordinateVector.fromValues(1,1,0);
		vertices[7] = CoordinateVector.fromValues(0,1,0);
		assertThrows(IllegalArgumentException.class, ()->new DistortedCell(vertices));
	}
	
	@Test
	public void testReferencePosition2D()
	{
		
		CoordinateVector [] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0,0);
		vertices[1] = CoordinateVector.fromValues(1,-0.5);
		vertices[2] = CoordinateVector.fromValues(14,1);
		vertices[3] = CoordinateVector.fromValues(-4.5, 6);
		DistortedCell cell = new DistortedCell(vertices);
		for(int i = 0; i < vertices.length; i++)
			assertEquals(i, cell.getPositionOfVertex(vertices[i]));
	}
	
	@Test
	public void testReferencePosition3D()
	{
		CoordinateVector [] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0,0,0);
		vertices[1] = CoordinateVector.fromValues(1,0,0);
		vertices[2] = CoordinateVector.fromValues(1,1,0);
		vertices[3] = CoordinateVector.fromValues(0,1,0);
		vertices[4] = CoordinateVector.fromValues(0,0,1);
		vertices[5] = CoordinateVector.fromValues(1,0,1);
		vertices[6] = CoordinateVector.fromValues(1,1,1);
		vertices[7] = CoordinateVector.fromValues(0,1,1);
		DistortedCell cell = new DistortedCell(vertices);
		for(int i = 0; i < vertices.length; i++)
			assertEquals(i, cell.getPositionOfVertex(vertices[i]));
	
	}
	
	@Test
	public void testGetters()
	{
	}
	
	@Test
	public void isInCell()
	{
	}
	
	@Test
	public void testNormal()
	{
	}
	
	@Test
	public void testCompare()
	{
	}
	
	@Test
	public void testEquals()
	{
	
	}
	
	@Test
	public void testHashCode()
	{
	
	}
	
	@Test
	public void testReferenceCell()
	{
	
	}
}
