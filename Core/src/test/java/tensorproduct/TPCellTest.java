package tensorproduct;

import basic.DoubleCompare;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.junit.jupiter.api.Test;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;


import static org.junit.jupiter.api.Assertions.*;

public class TPCellTest
{
	private TPCell createCell2D()
	{
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0,1), CoordinateVector.fromValues(6,
			1.0005), new IntCoordinates(1,1));
		return g.cells.get(0);
	}
	private TPCell createCell3D()
	{
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0,1,1676.2),
			CoordinateVector.fromValues(6,3,1889.4), new IntCoordinates(1,1,1));
		return g.cells.get(0);
	}
	@Test
	public void testGetters()
	{
		TPCell cell2d = createCell2D();
		assertEquals(4, cell2d.getFaces().size());
		assertEquals(2,cell2d.getDimension());
		assertEquals(CoordinateVector.fromValues(3,1.00025), cell2d.center());
		TPCell cell3d = createCell3D();
		assertEquals(6, cell3d.getFaces().size());
		assertEquals(3, cell3d.getDimension());
		assertTrue(CoordinateVector.fromValues(3,2, 1782.8).almostEqual(cell3d.center()));
	}
	@Test
	public void isInCell()
	{
		TPCell cell2d = createCell2D();
		assertTrue(cell2d.isInCell(CoordinateVector.fromValues(0,1)));
		assertTrue(cell2d.isInCell(CoordinateVector.fromValues(0,1.0002)));
		assertTrue(cell2d.isInCell(CoordinateVector.fromValues(0,1.0005)));
		assertTrue(cell2d.isInCell(CoordinateVector.fromValues(3,1)));
		assertTrue(cell2d.isInCell(CoordinateVector.fromValues(3,1.0002)));
		assertTrue(cell2d.isInCell(CoordinateVector.fromValues(3,1.0005)));
		assertTrue(cell2d.isInCell(CoordinateVector.fromValues(6,1)));
		assertTrue(cell2d.isInCell(CoordinateVector.fromValues(6,1.0002)));
		assertTrue(cell2d.isInCell(CoordinateVector.fromValues(6,1.0005)));
		assertFalse(cell2d.isInCell(CoordinateVector.fromValues(-1,1.0003)));
		assertFalse(cell2d.isInCell(CoordinateVector.fromValues(6.00001,1.0003)));
		assertFalse(cell2d.isInCell(CoordinateVector.fromValues(1,0)));
		assertFalse(cell2d.isInCell(CoordinateVector.fromValues(1,1.0005001)));
	}
	@Test
	public void testNormal()
	{
		CartesianGrid unrelated = new CartesianGrid(CoordinateVector.fromValues(5,6),
			CoordinateVector.fromValues(7,8), new IntCoordinates(1,1));
		TPCell cell2d = createCell2D();
		assertThrows(IllegalArgumentException.class,
			()->cell2d.getOuterNormal(unrelated.faces.get(0)));
		assertThrows(IllegalArgumentException.class,
			()->cell2d.getOuterNormal(unrelated.faces.get(1)));
		for(TPFace f:cell2d.getFaces())
		{
			CoordinateVector expectedNormal = new CoordinateVector(2);
			expectedNormal.set(f.otherCoordinate>1.0002?1:-1, f.flatDimension);
			assertEquals(expectedNormal, cell2d.getOuterNormal(f).value(f.center()));
		}
	}
	@Test
	public void testCompare()
	{
		TPCell cell1 = createCell2D();
		TPCell cell2 = createCell2D();
		assertEquals(0, cell1.compareTo(cell2));
	}
	@Test
	public void testGrid()
	{
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0,-1,1),
			CoordinateVector.fromValues(5,6,10), new IntCoordinates(5,7,10));
		for(TPCell cell: g.cells)
		{
			for (int i = 0; i < 3; i++)
			{
				assertTrue(DoubleCompare.isInteger(cell.getComponentCell(i).getStart()));
				assertTrue(DoubleCompare.isInteger(cell.getComponentCell(i).getEnd()));
			}
			assertEquals(6, cell.getFaces().size());
			for(TPFace f: cell.getFaces())
			{
				assertTrue((f.getCells().size() == 1 && f.isBoundaryFace()) || f.getCells().size() == 2);
			}
		}
		assertEquals(5*7*10, g.cells.size());
		assertEquals(6*7*10 + 5*8*10 + 5*7*11, g.faces.size());
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
	@Test
	public void testTransformationToReferenceCell()
	{
	
	}
	@Test
	public void testTransformationFromReferenceCell()
	{
	
	}
}