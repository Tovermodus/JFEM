package tensorproduct;

import basic.VectorFunction;
import linalg.CoordinateVector;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

public class TPCellTest
{
	private TPCell createCell2D()
	{
		List<Cell1D> cell1DS = List.of(new Cell1D(0,6), new Cell1D(1,1.0005));
		TPCell cell = new TPCell(cell1DS);
		cell.addFace(new TPFace(List.of(cell1DS.get(0)),1,1,true));
		cell.addFace(new TPFace(List.of(cell1DS.get(0)),1,1.0005,true));
		cell.addFace(new TPFace(List.of(cell1DS.get(1)),0,0,true));
		cell.addFace(new TPFace(List.of(cell1DS.get(1)),0,6,true));
		//cell.addFace(new TPFace(List.of(cell1DS.get(0)),1,1,true));
		return cell;
	}
	private TPCell createCell3D()
	{
		List<Cell1D> cell1DS = List.of(new Cell1D(0,6), new Cell1D(1,3), new Cell1D(1676.2, 1889.4));
		return new TPCell(cell1DS);
	}
	@Test
	public void testGetters()
	{
		TPCell cell2d = createCell2D();
		assertEquals(4, cell2d.getFaces().size());
		assertEquals(2,cell2d.getDimension());
		assertEquals(2,cell2d.getCell1Ds().size());
		assertEquals(CoordinateVector.fromValues(3,1.00025), cell2d.center());
		TPCell cell3d = createCell3D();
		assertEquals(0, cell3d.getFaces().size());
		assertEquals(3,cell3d.getDimension());
		assertEquals(3,cell3d.getCell1Ds().size());
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
		assertFalse(cell2d.isInCell(CoordinateVector.fromValues(1,1.0005000001)));
	}
	@Test
	public void testNormal()
	{
		TPCell cell2d = createCell2D();
		assertThrows(IllegalArgumentException.class,
			()->cell2d.getOuterNormal(new TPFace(List.of(cell2d.cell1Ds.get(0)),1,1,
				true)));
		assertThrows(IllegalArgumentException.class,
			()->cell2d.getOuterNormal(new TPFace(List.of(cell2d.cell1Ds.get(0)),0,1,
				true)));
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