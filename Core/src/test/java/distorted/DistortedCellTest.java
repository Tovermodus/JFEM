package distorted;

import distorted.geometry.CircleGrid;
import distorted.geometry.DistortedCell;
import linalg.CoordinateVector;
import org.junit.Test;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

public class DistortedCellTest
{
	private static DistortedCell createCell2D()
	{
		
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(1, -0.5);
		vertices[2] = CoordinateVector.fromValues(14, 1);
		vertices[3] = CoordinateVector.fromValues(-4.5, 6);
		final DistortedCell cell = new DistortedCell(vertices);
		CircleGrid.createFaces(new HashSet<>(), cell, List.of(cell), 100, vertices[0]);
		return cell;
	}
	
	private static DistortedCell createCell3D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0, 0, 3);
		vertices[1] = CoordinateVector.fromValues(1, 0, 3);
		vertices[2] = CoordinateVector.fromValues(1, 3, 3);
		vertices[3] = CoordinateVector.fromValues(0, 3, 3);
		vertices[4] = CoordinateVector.fromValues(0, 0, 4);
		vertices[5] = CoordinateVector.fromValues(2, 0, 4);
		vertices[6] = CoordinateVector.fromValues(2, 3, 4);
		vertices[7] = CoordinateVector.fromValues(0, 3, 4);
		final DistortedCell cell = new DistortedCell(vertices);
		CircleGrid.createFaces(new HashSet<>(), cell, List.of(cell), 100, vertices[0]);
		return cell;
	}
	
	@Test
	public void testOrdering2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 1);
		vertices[2] = CoordinateVector.fromValues(0, 1);
		vertices[3] = CoordinateVector.fromValues(1, 0);
		assertThrows(IllegalArgumentException.class, () -> new DistortedCell(vertices));
		
		vertices[1] = CoordinateVector.fromValues(0, 0);
		vertices[2] = CoordinateVector.fromValues(1, 0);
		vertices[3] = CoordinateVector.fromValues(1, 1);
		vertices[0] = CoordinateVector.fromValues(0, 1);
		assertThrows(IllegalArgumentException.class, () -> new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(0.5, 1.5);
		vertices[2] = CoordinateVector.fromValues(2, 2);
		vertices[3] = CoordinateVector.fromValues(0, 2);
		assertThrows(IllegalArgumentException.class, () -> new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(0, 2);
		vertices[2] = CoordinateVector.fromValues(2, 2);
		vertices[3] = CoordinateVector.fromValues(2, 0);
		assertThrows(IllegalArgumentException.class, () -> new DistortedCell(vertices));
	}
	
	@Test
	public void testOrdering3D()
	{
		
		final CoordinateVector[] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(-1, 0, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1, 0);
		vertices[3] = CoordinateVector.fromValues(0, 1, 0);
		vertices[4] = CoordinateVector.fromValues(0, 0, 1);
		vertices[5] = CoordinateVector.fromValues(1, 0, 1);
		vertices[6] = CoordinateVector.fromValues(1, 1, 1);
		vertices[7] = CoordinateVector.fromValues(0, 1, 1);
		assertThrows(IllegalArgumentException.class, () -> new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0, 2);
		vertices[2] = CoordinateVector.fromValues(1, 1, 0);
		vertices[3] = CoordinateVector.fromValues(0, 1, 0);
		vertices[4] = CoordinateVector.fromValues(0, 0, 1);
		vertices[5] = CoordinateVector.fromValues(1, 0, 1);
		vertices[6] = CoordinateVector.fromValues(1, 1, 1);
		vertices[7] = CoordinateVector.fromValues(0, 1, 1);
		assertThrows(IllegalArgumentException.class, () -> new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1, 0);
		vertices[3] = CoordinateVector.fromValues(0, 1, 0);
		vertices[4] = CoordinateVector.fromValues(0, 0, 1);
		vertices[5] = CoordinateVector.fromValues(0, 1, 1);
		vertices[6] = CoordinateVector.fromValues(1, 1, 1);
		vertices[7] = CoordinateVector.fromValues(1, 0, 1);
		assertThrows(IllegalArgumentException.class, () -> new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1, 0);
		vertices[3] = CoordinateVector.fromValues(0, 1, 0);
		vertices[4] = CoordinateVector.fromValues(0, 0, 1);
		vertices[5] = CoordinateVector.fromValues(1, 0, 1);
		vertices[6] = CoordinateVector.fromValues(1, 1, -1);
		vertices[7] = CoordinateVector.fromValues(0, 1, 1);
		assertThrows(IllegalArgumentException.class, () -> new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1, 0);
		vertices[3] = CoordinateVector.fromValues(0, 1, 0);
		vertices[4] = CoordinateVector.fromValues(1, 0, 1);
		vertices[5] = CoordinateVector.fromValues(1, 1, 1);
		vertices[6] = CoordinateVector.fromValues(0, 1, 1);
		vertices[7] = CoordinateVector.fromValues(0, 0, 1);
		assertThrows(IllegalArgumentException.class, () -> new DistortedCell(vertices));
		
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1, 1);
		vertices[3] = CoordinateVector.fromValues(0, 1, 1);
		vertices[4] = CoordinateVector.fromValues(0, 0, 1);
		vertices[5] = CoordinateVector.fromValues(1, 0, 1);
		vertices[6] = CoordinateVector.fromValues(1, 1, 0);
		vertices[7] = CoordinateVector.fromValues(0, 1, 0);
		assertThrows(IllegalArgumentException.class, () -> new DistortedCell(vertices));
	}
	
	@Test
	public void testReferencePosition2D()
	{
		
		final DistortedCell cell = createCell2D();
		assertEquals(CoordinateVector.fromValues(0, 0), cell.mapPositionToReferencePosition(0));
		assertEquals(CoordinateVector.fromValues(1, 0), cell.mapPositionToReferencePosition(1));
		assertEquals(CoordinateVector.fromValues(1, 1), cell.mapPositionToReferencePosition(2));
		assertEquals(CoordinateVector.fromValues(0, 1), cell.mapPositionToReferencePosition(3));
	}
	
	@Test
	public void testReferencePosition3D()
	{
		final DistortedCell cell = createCell3D();
		assertEquals(CoordinateVector.fromValues(0, 0, 0), cell.mapPositionToReferencePosition(0));
		assertEquals(CoordinateVector.fromValues(1, 0, 0), cell.mapPositionToReferencePosition(1));
		assertEquals(CoordinateVector.fromValues(1, 1, 0), cell.mapPositionToReferencePosition(2));
		assertEquals(CoordinateVector.fromValues(0, 1, 0), cell.mapPositionToReferencePosition(3));
		assertEquals(CoordinateVector.fromValues(0, 0, 1), cell.mapPositionToReferencePosition(4));
		assertEquals(CoordinateVector.fromValues(1, 0, 1), cell.mapPositionToReferencePosition(5));
		assertEquals(CoordinateVector.fromValues(1, 1, 1), cell.mapPositionToReferencePosition(6));
		assertEquals(CoordinateVector.fromValues(0, 1, 1), cell.mapPositionToReferencePosition(7));
	}
	
	@Test
	public void testVertexNumbers2D()
	{
		final DistortedCell cell = createCell2D();
		final List<CoordinateVector> vertices = cell.getVertices();
		for (int i = 0; i < vertices.size(); i++)
			assertEquals(i, cell.getPositionOfVertex(vertices.get(i)));
	}
	
	@Test
	public void testVertexNumbers3D()
	{
		
		final DistortedCell cell = createCell3D();
		final List<CoordinateVector> vertices = cell.getVertices();
		
		for (int i = 0; i < vertices.size(); i++)
			assertEquals(i, cell.getPositionOfVertex(vertices.get(i)));
	}
	
	@Test
	public void testGetters()
	{
		final DistortedCell cell2D = createCell2D();
		final DistortedCell cell3D = createCell3D();
		assertEquals(2, cell2D.getDimension());
		assertEquals(3, cell3D.getDimension());
	}
	
	@Test
	public void isInCell()
	{
		final DistortedCell cell2D = createCell2D();
		final DistortedCell cell3D = createCell3D();
		assertTrue(cell2D.isInCellPrecise(CoordinateVector.fromValues(0, 0)));
		assertTrue(cell2D.isInCellPrecise(CoordinateVector.fromValues(1, -0.5)));
		assertTrue(cell2D.isInCellPrecise(CoordinateVector.fromValues(14, 1)));
		assertTrue(cell2D.isInCellPrecise(CoordinateVector.fromValues(-4.5, 6)));
		assertTrue(cell2D.isInCellPrecise(CoordinateVector.fromValues(0.5, -0.25)));
		assertTrue(cell2D.isInCellPrecise(CoordinateVector.fromValues(7.5, 0.25)));
		assertTrue(cell2D.isInCellPrecise(CoordinateVector.fromValues(4.75, 2.5)));
		assertTrue(cell2D.isInCellPrecise(CoordinateVector.fromValues(-2.25, 3)));
		assertTrue(cell2D.isInCell(CoordinateVector.fromValues(1, 1)));
		
		assertFalse(cell2D.isInCell(CoordinateVector.fromValues(7.6, 0.25)));
		assertFalse(cell2D.isInCell(CoordinateVector.fromValues(1, -0.55)));
		assertFalse(cell2D.isInCell(CoordinateVector.fromValues(14, 0.9)));
		assertFalse(cell2D.isInCell(CoordinateVector.fromValues(0, 17)));
		assertFalse(cell2D.isInCell(CoordinateVector.fromValues(-4.6, 6)));
		assertFalse(cell2D.isInCell(CoordinateVector.fromValues(-10, 10)));
		
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(0, 0, 3)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(1, 0, 3)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(1, 3, 3)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(0, 3, 3)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(0, 0, 4)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(2, 0, 4)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(2, 3, 4)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(0, 3, 4)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(0, 1.5, 3)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(0.5, 0, 3)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(1.5, 3, 3.5)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(0, 3, 3.5)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(0, 0, 3.5)));
		assertTrue(cell3D.isInCellPrecise(CoordinateVector.fromValues(1, 1.5, 4)));
		assertTrue(cell3D.isInCell(CoordinateVector.fromValues(1, 1.5, 3.5)));
		
		assertFalse(cell3D.isInCellPrecise(CoordinateVector.fromValues(0, 0, 2.9)));
		assertFalse(cell3D.isInCellPrecise(CoordinateVector.fromValues(1.1, 0, 3)));
		assertFalse(cell3D.isInCellPrecise(CoordinateVector.fromValues(1, 3.1, 3)));
		assertFalse(cell3D.isInCellPrecise(CoordinateVector.fromValues(-0.1, 3, 3)));
		assertFalse(cell3D.isInCellPrecise(CoordinateVector.fromValues(0, 0, 4.1)));
		assertFalse(cell3D.isInCellPrecise(CoordinateVector.fromValues(2.1, 0, 4)));
		assertFalse(cell3D.isInCellPrecise(CoordinateVector.fromValues(2, 3.1, 4)));
		assertFalse(cell3D.isInCellPrecise(CoordinateVector.fromValues(-0.1, 3, 4)));
	}
	
	@Test
	public void testNormal()
	{
		final DistortedCell cell2D = createCell2D();
		final DistortedCell cell3D = createCell3D();
		final Set<CoordinateVector> normals2D = new HashSet<>(
			List.of(CoordinateVector.fromValues(-1, -2).normalize(),
			        CoordinateVector.fromValues(1.5, -13).normalize(),
			        CoordinateVector.fromValues(5, 18.5).normalize(),
			        CoordinateVector.fromValues(-6, -4.5).normalize()));
		final Set<CoordinateVector> normals3D = new HashSet<>(
			List.of(CoordinateVector.fromValues(0, 0, 1).normalize(),
			        CoordinateVector.fromValues(0, 0, -1).normalize(),
			        CoordinateVector.fromValues(0, -1, 0).normalize(),
			        CoordinateVector.fromValues(0, 1, 0).normalize(),
			        CoordinateVector.fromValues(-1, 0, 0).normalize(),
			        CoordinateVector.fromValues(1, 0, -1).normalize()));
		assertEquals(cell2D
			             .getFaces()
			             .stream()
			             .map(f -> cell2D.getOuterNormal(f).value(f.center()))
			             .collect(Collectors.toSet()), normals2D);
		assertEquals(cell3D
			             .getFaces()
			             .stream()
			             .map(f -> cell3D.getOuterNormal(f).value(f.center()))
			             .collect(Collectors.toSet()), normals3D);
	}
	
	@Test
	public void testCenter()
	{
		final DistortedCell cell2D = createCell2D();
		final DistortedCell cell3D = createCell3D();
		final DistortedCell cell2Dref = createCell2D().getReferenceCell();
		final DistortedCell cell3Dref = createCell3D().getReferenceCell();
		assertEquals(CoordinateVector.fromValues(0.5, 0.5), cell2Dref.center());
		assertEquals(CoordinateVector.fromValues(0.5, 0.5, 0.5), cell3Dref.center());
		assertEquals(CoordinateVector.fromValues(2.625, 1.625), cell2D.center());
		assertEquals(CoordinateVector.fromValues(0.75, 1.5, 3.5), cell3D.center());
	}
	
	@Test
	public void testEqualsCompareHashCode()
	{
		
		final DistortedCell reference = createCell2D();
		CoordinateVector[] vertices = new CoordinateVector[reference.getVertices().size()];
		DistortedCell other;
		for (int i = 0; i < reference.getVertices().size(); i++)
		{
			vertices = reference.getVertices().toArray(vertices).clone();
			vertices[i].addInPlace(CoordinateVector.fromValues(0, 1e-5));
			other = new DistortedCell(vertices);
			assertNotEquals(reference, other);
			assertTrue(reference.compareTo(other) != 0);
			assertTrue(other.compareTo(reference) != 0);
		}
		for (int i = 0; i < reference.getVertices().size(); i++)
		{
			vertices = reference.getVertices().toArray(vertices);
			vertices[i].addInPlace(CoordinateVector.fromValues(0, 1e-13));
			other = new DistortedCell(vertices);
			assertEquals(reference.hashCode(), other.hashCode());
			assertEquals(reference.hashCode(), other.hashCode());
			assertEquals(0, reference.compareTo(other));
			assertEquals(0, other.compareTo(reference));
		}
	}
}
