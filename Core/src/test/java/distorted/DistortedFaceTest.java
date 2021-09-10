package distorted;

import distorted.geometry.CircleGrid;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.AffineTransformation;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.junit.Test;

import java.util.HashSet;
import java.util.List;
import java.util.Optional;

import static org.junit.Assert.*;

public class DistortedFaceTest
{
	private static DistortedFace getCentralFace2D(final int side)
	{
		final HashSet<DistortedFace> faces = createFaces2D(side);
		final Optional<DistortedFace> centralFace =
			faces.stream().filter(f -> f.getCells().size() == 2).findAny();
		if (centralFace.isPresent())
			return centralFace.get();
		else throw new IllegalArgumentException("Wrong side");
	}
	
	private static HashSet<DistortedFace> createFaces2D(final int side)
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(1, -0.5);
		vertices[2] = CoordinateVector.fromValues(14, 1);
		vertices[3] = CoordinateVector.fromValues(-4.5, 6);
		final DistortedCell cell = new DistortedCell(vertices);
		final CoordinateVector[] otherVertices = new CoordinateVector[4];
		if (side == 0)
		{
			otherVertices[0] = CoordinateVector.fromValues(-1, -1);
			otherVertices[1] = CoordinateVector.fromValues(1, -1);
			otherVertices[2] = CoordinateVector.fromValues(1, -0.5);
			otherVertices[3] = CoordinateVector.fromValues(0, 0);
		}
		if (side == 1)
		{
			otherVertices[0] = CoordinateVector.fromValues(1, -1);
			otherVertices[1] = CoordinateVector.fromValues(5, -1);
			otherVertices[2] = CoordinateVector.fromValues(14, 1);
			otherVertices[3] = CoordinateVector.fromValues(1, -0.5);
		}
		if (side == 2)
		{
			otherVertices[0] = CoordinateVector.fromValues(-4.5, 6);
			otherVertices[1] = CoordinateVector.fromValues(14, 1);
			otherVertices[2] = CoordinateVector.fromValues(14, 5);
			otherVertices[3] = CoordinateVector.fromValues(-4.5, 10);
		}
		if (side == 3)
		{
			otherVertices[0] = CoordinateVector.fromValues(-1, 0);
			otherVertices[1] = CoordinateVector.fromValues(0, 0);
			otherVertices[2] = CoordinateVector.fromValues(-4.5, 6);
			otherVertices[3] = CoordinateVector.fromValues(-7, 8);
		}
		final DistortedCell otherCell = new DistortedCell(otherVertices);
		final HashSet<DistortedFace> faces = new HashSet<>();
		CircleGrid.createFaces(faces, cell, List.of(cell, otherCell),
		                       vertices[0].sub(vertices[1]).euclidianNorm() / 2,
		                       vertices[0].add(vertices[1]).mul(0.5));
		CircleGrid.createFaces(faces, otherCell, List.of(cell, otherCell),
		                       vertices[0].sub(vertices[1]).euclidianNorm() / 2,
		                       vertices[0].add(vertices[1]).mul(0.5));
		return faces;
	}
	
	@Test
	public void testRotation2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1);
		vertices[3] = CoordinateVector.fromValues(0, 1);
		final DistortedCell cell = new DistortedCell(vertices);
		for (int i = 0; i < vertices.length; i++)
		{
			final DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[i],
			                                                                    vertices[(i + 1) %
				                                                                    vertices.length]},
			                                             false);
			final AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(
				cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0)), vertices[i]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0)), vertices[(i + 1) % vertices.length]);
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(vertices[(i + 2) % vertices.length]));
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(vertices[(i + 3) % vertices.length]));
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(vertices[(i) % vertices.length]));
		}
		for (int i = 0; i < vertices.length; i++)
		{
			final DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[i],
			                                                                    vertices[(i + 3) %
				                                                                    vertices.length]},
			                                             false);
			final AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(
				cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0)), vertices[i]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0)), vertices[(i + 3) % vertices.length]);
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(vertices[(i) % vertices.length]));
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(vertices[(i + 1) % vertices.length]));
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(vertices[(i + 2) % vertices.length]));
		}
	}
	
	@Test
	public void testDistortedRotation2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(2, 1);
		vertices[2] = CoordinateVector.fromValues(4, 4);
		vertices[3] = CoordinateVector.fromValues(-0.5, 2);
		final DistortedCell cell = new DistortedCell(vertices);
		final CoordinateVector[] referenceVertices = new CoordinateVector[4];
		referenceVertices[0] = CoordinateVector.fromValues(0, 0);
		referenceVertices[1] = CoordinateVector.fromValues(1, 0);
		referenceVertices[2] = CoordinateVector.fromValues(1, 1);
		referenceVertices[3] = CoordinateVector.fromValues(0, 1);
		
		for (int i = 0; i < vertices.length; i++)
		{
			final DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[i],
			                                                                    vertices[(i + 1) %
				                                                                    vertices.length]},
			                                             false);
			final AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(
				cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0)), referenceVertices[i]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0)),
			             referenceVertices[(i + 1) % vertices.length]);
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(referenceVertices[(i + 2) % vertices.length]));
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(referenceVertices[(i + 3) % vertices.length]));
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(referenceVertices[(i) % vertices.length]));
		}
		for (int i = 0; i < vertices.length; i++)
		{
			final DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[i],
			                                                                    vertices[(i + 3) %
				                                                                    vertices.length]},
			                                             false);
			final AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(
				cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0)), referenceVertices[i]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0)),
			             referenceVertices[(i + 3) % vertices.length]);
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(referenceVertices[(i) % vertices.length]));
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(referenceVertices[(i + 1) % vertices.length]));
			assertFalse(at
				            .apply(CoordinateVector.fromValues(1, 0))
				            .almostEqual(referenceVertices[(i + 2) % vertices.length]));
		}
	}
	
	@Test
	public void testRotation3D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1, 0);
		vertices[3] = CoordinateVector.fromValues(0, 1, 0);
		vertices[4] = CoordinateVector.fromValues(0, 0, 1);
		vertices[5] = CoordinateVector.fromValues(1, 0, 1);
		vertices[6] = CoordinateVector.fromValues(1, 1, 1);
		vertices[7] = CoordinateVector.fromValues(0, 1, 1);
		final DistortedCell cell = new DistortedCell(vertices);
		final IntCoordinates[] faces = new IntCoordinates[48];
		faces[0] = new IntCoordinates(0, 3, 2, 1);
		faces[1] = new IntCoordinates(3, 2, 1, 0);
		faces[2] = new IntCoordinates(2, 1, 0, 3);
		faces[3] = new IntCoordinates(1, 0, 3, 2);
		faces[4] = new IntCoordinates(0, 1, 2, 3);
		faces[5] = new IntCoordinates(1, 2, 3, 0);
		faces[6] = new IntCoordinates(2, 3, 0, 1);
		faces[7] = new IntCoordinates(3, 0, 1, 2);
		
		faces[8] = new IntCoordinates(0, 1, 5, 4);
		faces[9] = new IntCoordinates(1, 5, 4, 0);
		faces[10] = new IntCoordinates(5, 4, 0, 1);
		faces[11] = new IntCoordinates(4, 0, 1, 5);
		faces[12] = new IntCoordinates(0, 4, 5, 1);
		faces[13] = new IntCoordinates(4, 5, 1, 0);
		faces[14] = new IntCoordinates(5, 1, 0, 4);
		faces[15] = new IntCoordinates(1, 0, 4, 5);
		
		faces[16] = new IntCoordinates(0, 4, 7, 3);
		faces[17] = new IntCoordinates(4, 7, 3, 0);
		faces[18] = new IntCoordinates(7, 3, 0, 4);
		faces[19] = new IntCoordinates(3, 0, 4, 7);
		faces[20] = new IntCoordinates(0, 3, 7, 4);
		faces[21] = new IntCoordinates(3, 7, 4, 0);
		faces[22] = new IntCoordinates(7, 4, 0, 3);
		faces[23] = new IntCoordinates(4, 0, 3, 7);
		
		faces[24] = new IntCoordinates(6, 7, 4, 5);
		faces[25] = new IntCoordinates(7, 4, 5, 6);
		faces[26] = new IntCoordinates(4, 5, 6, 7);
		faces[27] = new IntCoordinates(5, 6, 7, 4);
		faces[28] = new IntCoordinates(6, 5, 4, 7);
		faces[29] = new IntCoordinates(5, 4, 7, 6);
		faces[30] = new IntCoordinates(4, 7, 6, 5);
		faces[31] = new IntCoordinates(7, 6, 5, 4);
		
		faces[32] = new IntCoordinates(6, 2, 3, 7);
		faces[33] = new IntCoordinates(2, 3, 7, 6);
		faces[34] = new IntCoordinates(3, 7, 6, 2);
		faces[35] = new IntCoordinates(7, 6, 2, 3);
		faces[36] = new IntCoordinates(6, 7, 3, 2);
		faces[37] = new IntCoordinates(7, 3, 2, 6);
		faces[38] = new IntCoordinates(3, 2, 6, 7);
		faces[39] = new IntCoordinates(2, 6, 7, 3);
		
		faces[40] = new IntCoordinates(6, 5, 1, 2);
		faces[41] = new IntCoordinates(5, 1, 2, 6);
		faces[42] = new IntCoordinates(1, 2, 6, 5);
		faces[43] = new IntCoordinates(2, 6, 5, 1);
		faces[44] = new IntCoordinates(6, 2, 1, 5);
		faces[45] = new IntCoordinates(2, 1, 5, 6);
		faces[46] = new IntCoordinates(1, 5, 6, 2);
		faces[47] = new IntCoordinates(5, 6, 2, 1);
		
		for (final IntCoordinates facec : faces)
		{
			final DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[facec.get(0)],
			                                                                    vertices[facec.get(1)],
			                                                                    vertices[facec.get(2)],
			                                                                    vertices[facec.get(3)]},
			                                             false);
			final AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(
				cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0, 0)), vertices[facec.get(0)]);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 1, 0)), vertices[facec.get(1)]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 1, 0)), vertices[facec.get(2)]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0, 0)), vertices[facec.get(3)]);
		}
	}
	
	@Test
	public void testDistortedRotation3D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1, 0);
		vertices[3] = CoordinateVector.fromValues(0, 1, 0);
		vertices[4] = CoordinateVector.fromValues(0, 0, 1);
		vertices[5] = CoordinateVector.fromValues(1, 0, 1);
		vertices[6] = CoordinateVector.fromValues(1, 1, 1);
		vertices[7] = CoordinateVector.fromValues(0, 1, 1);
		final CoordinateVector[] referenceVertices = new CoordinateVector[8];
		referenceVertices[0] = CoordinateVector.fromValues(0, 0, 0);
		referenceVertices[1] = CoordinateVector.fromValues(1, 0, 0);
		referenceVertices[2] = CoordinateVector.fromValues(1, 1, 0);
		referenceVertices[3] = CoordinateVector.fromValues(0, 1, 0);
		referenceVertices[4] = CoordinateVector.fromValues(0, 0, 1);
		referenceVertices[5] = CoordinateVector.fromValues(1, 0, 1);
		referenceVertices[6] = CoordinateVector.fromValues(1, 1, 1);
		referenceVertices[7] = CoordinateVector.fromValues(0, 1, 1);
		final DistortedCell cell = new DistortedCell(vertices);
		final IntCoordinates[] faces = new IntCoordinates[48];
		faces[0] = new IntCoordinates(0, 3, 2, 1);
		faces[1] = new IntCoordinates(3, 2, 1, 0);
		faces[2] = new IntCoordinates(2, 1, 0, 3);
		faces[3] = new IntCoordinates(1, 0, 3, 2);
		faces[4] = new IntCoordinates(0, 1, 2, 3);
		faces[5] = new IntCoordinates(1, 2, 3, 0);
		faces[6] = new IntCoordinates(2, 3, 0, 1);
		faces[7] = new IntCoordinates(3, 0, 1, 2);
		
		faces[8] = new IntCoordinates(0, 1, 5, 4);
		faces[9] = new IntCoordinates(1, 5, 4, 0);
		faces[10] = new IntCoordinates(5, 4, 0, 1);
		faces[11] = new IntCoordinates(4, 0, 1, 5);
		faces[12] = new IntCoordinates(0, 4, 5, 1);
		faces[13] = new IntCoordinates(4, 5, 1, 0);
		faces[14] = new IntCoordinates(5, 1, 0, 4);
		faces[15] = new IntCoordinates(1, 0, 4, 5);
		
		faces[16] = new IntCoordinates(0, 4, 7, 3);
		faces[17] = new IntCoordinates(4, 7, 3, 0);
		faces[18] = new IntCoordinates(7, 3, 0, 4);
		faces[19] = new IntCoordinates(3, 0, 4, 7);
		faces[20] = new IntCoordinates(0, 3, 7, 4);
		faces[21] = new IntCoordinates(3, 7, 4, 0);
		faces[22] = new IntCoordinates(7, 4, 0, 3);
		faces[23] = new IntCoordinates(4, 0, 3, 7);
		
		faces[24] = new IntCoordinates(6, 7, 4, 5);
		faces[25] = new IntCoordinates(7, 4, 5, 6);
		faces[26] = new IntCoordinates(4, 5, 6, 7);
		faces[27] = new IntCoordinates(5, 6, 7, 4);
		faces[28] = new IntCoordinates(6, 5, 4, 7);
		faces[29] = new IntCoordinates(5, 4, 7, 6);
		faces[30] = new IntCoordinates(4, 7, 6, 5);
		faces[31] = new IntCoordinates(7, 6, 5, 4);
		
		faces[32] = new IntCoordinates(6, 2, 3, 7);
		faces[33] = new IntCoordinates(2, 3, 7, 6);
		faces[34] = new IntCoordinates(3, 7, 6, 2);
		faces[35] = new IntCoordinates(7, 6, 2, 3);
		faces[36] = new IntCoordinates(6, 7, 3, 2);
		faces[37] = new IntCoordinates(7, 3, 2, 6);
		faces[38] = new IntCoordinates(3, 2, 6, 7);
		faces[39] = new IntCoordinates(2, 6, 7, 3);
		
		faces[40] = new IntCoordinates(6, 5, 1, 2);
		faces[41] = new IntCoordinates(5, 1, 2, 6);
		faces[42] = new IntCoordinates(1, 2, 6, 5);
		faces[43] = new IntCoordinates(2, 6, 5, 1);
		faces[44] = new IntCoordinates(6, 2, 1, 5);
		faces[45] = new IntCoordinates(2, 1, 5, 6);
		faces[46] = new IntCoordinates(1, 5, 6, 2);
		faces[47] = new IntCoordinates(5, 6, 2, 1);
		
		for (final IntCoordinates facec : faces)
		{
			final DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[facec.get(0)],
			                                                                    vertices[facec.get(1)],
			                                                                    vertices[facec.get(2)],
			                                                                    vertices[facec.get(3)]},
			                                             false);
			final AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(
				cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0, 0)), referenceVertices[facec.get(0)]);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 1, 0)), referenceVertices[facec.get(1)]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 1, 0)), referenceVertices[facec.get(2)]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0, 0)), referenceVertices[facec.get(3)]);
		}
	}
	
	@Test
	public void testIsOnCell()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(1, -0.5);
		vertices[2] = CoordinateVector.fromValues(14, 1);
		vertices[3] = CoordinateVector.fromValues(-4.5, 6);
		final DistortedCell on = new DistortedCell(vertices);
		vertices[0] = CoordinateVector.fromValues(0.1, 0);
		vertices[1] = CoordinateVector.fromValues(1.1, -0.5);
		vertices[2] = CoordinateVector.fromValues(14.1, 1);
		vertices[3] = CoordinateVector.fromValues(-4.6, 6);
		final DistortedCell off = new DistortedCell(vertices);
		for (int i = 0; i < 4; i++)
		{
			final DistortedFace f = getCentralFace2D(i);
			assertTrue(f.isOnCell(on));
			assertFalse(f.isOnCell(off));
		}
	}
	
	@Test
	public void testCenter()
	{
		for (int i = 0; i < 4; i++)
		{
			for (final DistortedFace f : createFaces2D(i))
			{
				if (f.getNormalUpstreamCell() != null)
					assertTrue(f.getNormalUpstreamCell().isInCellPrecise(f.center()));
				if (f.getNormalDownstreamCell() != null)
					assertTrue(f.getNormalDownstreamCell().isInCellPrecise(f.center()));
				assertTrue(f.isOnFace(f.center()));
			}
		}
	}
	
	@Test
	public void testGetDownStreamCell()
	{
		for (int i = 0; i < 4; i++)
		{
			final DistortedFace f = getCentralFace2D(i);
			assertNotEquals(CoordinateVector.fromValues(2.625, 1.625),
			                f.getNormalDownstreamCell().center());
		}
		for (int i = 0; i < 4; i++)
		{
			for (final DistortedFace f : createFaces2D(i))
			{
				if (f.getCells().size() != 2)
					assertNull(f.getNormalDownstreamCell());
				else
					assertNotNull(f.getNormalDownstreamCell());
			}
		}
	}
	
	@Test
	public void testGetUpStreamCell()
	{
		for (int i = 0; i < 4; i++)
		{
			final DistortedFace f = getCentralFace2D(i);
			assertEquals(CoordinateVector.fromValues(2.625, 1.625), f.getNormalUpstreamCell().center());
		}
		for (int i = 0; i < 4; i++)
		{
			for (final DistortedFace f : createFaces2D(i))
			{
				assertNotNull(f.getNormalUpstreamCell());
			}
		}
	}
	
	@Test
	public void testEqualsCompareHashCode()
	{
		for (int side = 0; side < 4; side++)
		{
			final DistortedFace reference = getCentralFace2D(side);
			CoordinateVector[] vertices = reference.getVertices().toArray(new CoordinateVector[2]);
			DistortedFace other = new DistortedFace(vertices, !reference.isBoundaryFace());
			assertNotEquals(reference.hashCode(), other.hashCode());
			assertNotEquals(reference, other);
			assertTrue(reference.compareTo(other) != 0);
			assertTrue(other.compareTo(reference) != 0);
			for (int i = 0; i < reference.getVertices().size(); i++)
			{
				vertices = reference.getVertices().toArray(vertices).clone();
				vertices[i].addInPlace(CoordinateVector.fromValues(0, 1e-5));
				other = new DistortedFace(vertices, reference.isBoundaryFace());
				assertNotEquals(reference, other);
				assertTrue(reference.compareTo(other) != 0);
				assertTrue(other.compareTo(reference) != 0);
			}
			for (int i = 0; i < reference.getVertices().size(); i++)
			{
				vertices = reference.getVertices().toArray(vertices);
				vertices[i].addInPlace(CoordinateVector.fromValues(0, 1e-13));
				other = new DistortedFace(vertices, reference.isBoundaryFace());
				assertEquals(reference.hashCode(), other.hashCode());
				assertEquals(reference.hashCode(), other.hashCode());
				assertEquals(0, reference.compareTo(other));
				assertEquals(0, other.compareTo(reference));
			}
		}
	}
	
	@Test
	public void testIsVertex()
	{
		final DistortedCell cell = getCentralFace2D(0).getNormalUpstreamCell();
		for (int i = 0; i < 4; i++)
		{
			for (final DistortedFace f : cell.getFaces())
			{
				assertEquals(2, cell.getVertices().stream().filter(f::isVertex).count());
			}
		}
	}
}
