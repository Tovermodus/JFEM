package distorted;

import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.AffineTransformation;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;

public class DistortedFaceTransformationTest
{
	@Test
	public void testRotation2D()
	{
		CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1);
		vertices[3] = CoordinateVector.fromValues(0, 1);
		DistortedCell cell = new DistortedCell(vertices);
		for (int i = 0; i < vertices.length; i++)
		{
			DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[i],
				vertices[(i + 1) % vertices.length]}, false);
			AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0)), vertices[i]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0)), vertices[(i + 1) % vertices.length]);
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(vertices[(i + 2) % vertices.length]));
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(vertices[(i + 3) % vertices.length]));
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(vertices[(i) % vertices.length]));
			
		}
		for (int i = 0; i < vertices.length; i++)
		{
			DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[i],
				vertices[(i + 3) % vertices.length]}, false);
			AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0)), vertices[i]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0)), vertices[(i + 3) % vertices.length]);
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(vertices[(i) % vertices.length]));
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(vertices[(i + 1) % vertices.length]));
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(vertices[(i + 2) % vertices.length]));
			
		}
	}
	
	@Test
	public void testDistortedRotation2D()
	{
		CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(2, 1);
		vertices[2] = CoordinateVector.fromValues(4, 4);
		vertices[3] = CoordinateVector.fromValues(-0.5, 2);
		DistortedCell cell = new DistortedCell(vertices);
		CoordinateVector[] referenceVertices = new CoordinateVector[4];
		referenceVertices[0] = CoordinateVector.fromValues(0, 0);
		referenceVertices[1] = CoordinateVector.fromValues(1, 0);
		referenceVertices[2] = CoordinateVector.fromValues(1, 1);
		referenceVertices[3] = CoordinateVector.fromValues(0, 1);
		
		for (int i = 0; i < vertices.length; i++)
		{
			DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[i],
				vertices[(i + 1) % vertices.length]}, false);
			AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0)), referenceVertices[i]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0)), referenceVertices[(i + 1) % vertices.length]);
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(referenceVertices[(i + 2) % vertices.length]));
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(referenceVertices[(i + 3) % vertices.length]));
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(referenceVertices[(i) % vertices.length]));
			
		}
		for (int i = 0; i < vertices.length; i++)
		{
			DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[i],
				vertices[(i + 3) % vertices.length]}, false);
			AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0)), referenceVertices[i]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0)), referenceVertices[(i + 3) % vertices.length]);
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(referenceVertices[(i) % vertices.length]));
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(referenceVertices[(i + 1) % vertices.length]));
			assertFalse(at.apply(CoordinateVector.fromValues(1, 0)).almostEqual(referenceVertices[(i + 2) % vertices.length]));
			
		}
		
	}
	@Test
	public void testRotation3D()
	{
		CoordinateVector[] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1, 0);
		vertices[3] = CoordinateVector.fromValues(0, 1, 0);
		vertices[4] = CoordinateVector.fromValues(0, 0, 1);
		vertices[5] = CoordinateVector.fromValues(1, 0, 1);
		vertices[6] = CoordinateVector.fromValues(1, 1, 1);
		vertices[7] = CoordinateVector.fromValues(0, 1, 1);
		DistortedCell cell = new DistortedCell(vertices);
		IntCoordinates[] faces = new IntCoordinates[48];
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
		
		for (IntCoordinates facec : faces)
		{
			DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[facec.get(0)],
				vertices[facec.get(1)],
				vertices[facec.get(2)],
				vertices[facec.get(3)]}, false);
			AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0, 0)), vertices[facec.get(0)]);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 1, 0)), vertices[facec.get(1)]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 1, 0)), vertices[facec.get(2)]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0, 0)), vertices[facec.get(3)]);
			
		}
	}
	@Test
	public void testDistortedRotation3D()
	{
		CoordinateVector[] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1, 0);
		vertices[3] = CoordinateVector.fromValues(0, 1, 0);
		vertices[4] = CoordinateVector.fromValues(0, 0, 1);
		vertices[5] = CoordinateVector.fromValues(1, 0, 1);
		vertices[6] = CoordinateVector.fromValues(1, 1, 1);
		vertices[7] = CoordinateVector.fromValues(0, 1, 1);
		CoordinateVector[] referenceVertices = new CoordinateVector[8];
		referenceVertices[0] = CoordinateVector.fromValues(0, 0, 0);
		referenceVertices[1] = CoordinateVector.fromValues(1, 0, 0);
		referenceVertices[2] = CoordinateVector.fromValues(1, 1, 0);
		referenceVertices[3] = CoordinateVector.fromValues(0, 1, 0);
		referenceVertices[4] = CoordinateVector.fromValues(0, 0, 1);
		referenceVertices[5] = CoordinateVector.fromValues(1, 0, 1);
		referenceVertices[6] = CoordinateVector.fromValues(1, 1, 1);
		referenceVertices[7] = CoordinateVector.fromValues(0, 1, 1);
		DistortedCell cell = new DistortedCell(vertices);
		IntCoordinates[] faces = new IntCoordinates[48];
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
		
		for (IntCoordinates facec : faces)
		{
			DistortedFace face = new DistortedFace(new CoordinateVector[]{vertices[facec.get(0)],
				vertices[facec.get(1)],
				vertices[facec.get(2)],
				vertices[facec.get(3)]}, false);
			AffineTransformation at = face.getTransformationFromReferenceFaceToFaceOfReferenceCell(cell);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 0, 0)), referenceVertices[facec.get(0)]);
			assertEquals(at.apply(CoordinateVector.fromValues(0, 1, 0)), referenceVertices[facec.get(1)]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 1, 0)), referenceVertices[facec.get(2)]);
			assertEquals(at.apply(CoordinateVector.fromValues(1, 0, 0)), referenceVertices[facec.get(3)]);
			
		}
	}
}
