package distorted;

import distorted.geometry.DistortedCell;
import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class DistortedMappingTest
{
	private static DistortedCell getNonDistortedCell2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(1, 0);
		vertices[2] = CoordinateVector.fromValues(1, 1);
		vertices[3] = CoordinateVector.fromValues(0, 1);
		return new DistortedCell(vertices);
	}
	
	private static DistortedCell getNonDistortedCell3D()
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
		return new DistortedCell(vertices);
	}
	
	private static DistortedCell getScaled2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(4, 0);
		vertices[2] = CoordinateVector.fromValues(4, 2);
		vertices[3] = CoordinateVector.fromValues(0, 2);
		return new DistortedCell(vertices);
	}
	
	private static DistortedCell getScaled3D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(4, 0, 0);
		vertices[2] = CoordinateVector.fromValues(4, 2, 0);
		vertices[3] = CoordinateVector.fromValues(0, 2, 0);
		vertices[4] = CoordinateVector.fromValues(0, 0, 8);
		vertices[5] = CoordinateVector.fromValues(4, 0, 8);
		vertices[6] = CoordinateVector.fromValues(4, 2, 8);
		vertices[7] = CoordinateVector.fromValues(0, 2, 8);
		return new DistortedCell(vertices);
	}
	
	private static DistortedCell getRotated2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(1, -1);
		vertices[2] = CoordinateVector.fromValues(2, 0);
		vertices[3] = CoordinateVector.fromValues(1, 1);
		return new DistortedCell(vertices);
	}
	
	private static DistortedCell getRotated3D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0, 0, 0);
		vertices[1] = CoordinateVector.fromValues(1, -1, 0);
		vertices[2] = CoordinateVector.fromValues(2, 0, 0);
		vertices[3] = CoordinateVector.fromValues(1, 1, 0);
		vertices[4] = CoordinateVector.fromValues(0, 0, 1);
		vertices[5] = CoordinateVector.fromValues(1, -1, 1);
		vertices[6] = CoordinateVector.fromValues(2, 0, 1);
		vertices[7] = CoordinateVector.fromValues(1, 1, 1);
		return new DistortedCell(vertices);
	}
	
	private static DistortedCell getDistortedCell2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0, 0);
		vertices[1] = CoordinateVector.fromValues(2, 0);
		vertices[2] = CoordinateVector.fromValues(0.8, 23.1);
		vertices[3] = CoordinateVector.fromValues(0, 1);
		return new DistortedCell(vertices);
	}
	
	private static DistortedCell getAnotherDistortedCell2D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(-2, -2);
		vertices[1] = CoordinateVector.fromValues(0, -4 / Math.sqrt(2));
		vertices[2] = CoordinateVector.fromValues(0, -1 - 1. / Math.sqrt(2));
		vertices[3] = CoordinateVector.fromValues(-1.5, -1.5);
		return new DistortedCell(vertices);
	}
	
	private static DistortedCell getDistortedCell3D()
	{
		final CoordinateVector[] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(5, 5, 5);
		vertices[1] = CoordinateVector.fromValues(10, 0, 5);
		vertices[2] = CoordinateVector.fromValues(15, 5, 5);
		vertices[3] = CoordinateVector.fromValues(10, 10, 5);
		vertices[4] = CoordinateVector.fromValues(4, 6, 11);
		vertices[5] = CoordinateVector.fromValues(9, 1, 11);
		vertices[6] = CoordinateVector.fromValues(14, 6, 11);
		vertices[7] = CoordinateVector.fromValues(9, 11, 11);
		return new DistortedCell(vertices);
	}
	
	@Test
	public void testMappingNonDistorted3D()
	{
		final DistortedCell cell = getNonDistortedCell3D();
		for (int i = 0; i <= 10; i++)
		{
			for (int j = 0; j <= 10; j++)
			{
				for (int k = 0; k <= 10; k++)
				{
					final CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1,
					                                                       k * 0.1);
					assertTrue(c.almostEqual(cell.transformPreciseToReferenceCell(c), 1e-8));
					assertTrue(c.almostEqual(cell.transformFromReferenceCell(c), 1e-8));
					assertEquals(cell.transformationGradientFromReferenceCell(c),
					             CoordinateDenseMatrix.fromValues(3, 3, 1, 0, 0, 0, 1, 0, 0, 0, 1));
					assertEquals(cell.transformationGradientToReferenceCell(c),
					             CoordinateDenseMatrix.fromValues(3, 3, 1, 0, 0, 0, 1, 0, 0, 0, 1));
					assertTrue(c.almostEqual(cell.transformPreciseToReferenceCell(
						cell.transformFromReferenceCell(c)), 1e-8));
				}
			}
		}
	}
	
	@Test
	public void testMappingNonDistorted2D()
	{
		final DistortedCell cell = getNonDistortedCell2D();
		for (int i = 0; i <= 10; i++)
		{
			for (int j = 0; j <= 10; j++)
			{
				final CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1);
				assertTrue(c.almostEqual(cell.transformPreciseToReferenceCell(c), 1e-8));
				assertTrue(c.almostEqual(cell.transformFromReferenceCell(c), 1e-8));
				assertEquals(cell.transformationGradientFromReferenceCell(c),
				             CoordinateDenseMatrix.fromValues(2, 2, 1, 0, 0, 1));
				assertEquals(cell.transformationGradientToReferenceCell(c),
				             CoordinateDenseMatrix.fromValues(2, 2, 1, 0, 0, 1));
				assertTrue(c.almostEqual(
					cell.transformPreciseToReferenceCell(cell.transformFromReferenceCell(c)),
					1e-8));
			}
		}
	}
	
	@Test
	public void testMappingScaled2D()
	{
		final DistortedCell cell = getScaled2D();
		final CoordinateMatrix scale = CoordinateDenseMatrix.fromValues(2, 2, 4, 0, 0, 2);
		final CoordinateMatrix scaleInv = CoordinateDenseMatrix.fromValues(2, 2, 0.25, 0, 0, 0.5);
		for (int i = 0; i <= 10; i++)
		{
			for (int j = 0; j <= 10; j++)
			{
				final CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1);
				assertEquals(cell.transformationGradientFromReferenceCell(c), scale);
				assertEquals(cell.transformationGradientToReferenceCell(c), scaleInv);
				assertTrue(c.almostEqual(cell.transformPreciseToReferenceCell(scale.mvMul(c)), 1e-8));
				assertTrue(c.almostEqual(cell.transformFromReferenceCell(scaleInv.mvMul(c)), 1e-8));
				assertTrue(c.almostEqual(
					cell.transformPreciseToReferenceCell(cell.transformFromReferenceCell(c)),
					1e-8));
			}
		}
	}
	
	@Test
	public void testMappingScaled3D()
	{
		final DistortedCell cell = getScaled3D();
		final CoordinateMatrix scale = CoordinateDenseMatrix.fromValues(3, 3, 4, 0, 0, 0, 2, 0, 0, 0, 8);
		final CoordinateMatrix scaleInv = CoordinateDenseMatrix.fromValues(3, 3, 0.25, 0, 0, 0, 0.5, 0, 0, 0,
		                                                                   0.125);
		for (int i = 0; i <= 10; i++)
		{
			for (int j = 0; j <= 10; j++)
			{
				for (int k = 0; k <= 10; k++)
				{
					final CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1,
					                                                       k * 0.1);
					assertTrue(c.almostEqual(cell.transformPreciseToReferenceCell(scale.mvMul(c)),
					                         1e-8));
					assertTrue(c.almostEqual(cell.transformFromReferenceCell(scaleInv.mvMul(c)),
					                         1e-8));
					assertEquals(cell.transformationGradientFromReferenceCell(c), scale);
					assertEquals(cell.transformationGradientToReferenceCell(c), scaleInv);
					assertTrue(c.almostEqual(cell.transformPreciseToReferenceCell(
						cell.transformFromReferenceCell(c)), 1e-8));
				}
			}
		}
	}
	
	@Test
	public void testMappingRotated2D()
	{
		final DistortedCell cell = getRotated2D();
		final CoordinateDenseMatrix rotation = CoordinateDenseMatrix.fromValues(2, 2, 1, 1, -1, 1);
		final CoordinateMatrix rotationInv = rotation.inverse();
		for (int i = 0; i <= 10; i++)
		{
			for (int j = 0; j <= 10; j++)
			{
				final CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1);
				assertEquals(cell.transformationGradientFromReferenceCell(c), rotation);
				assertEquals(cell.transformationGradientToReferenceCell(c), rotationInv);
				assertTrue(c.almostEqual(
					cell.transformPreciseToReferenceCell(cell.transformFromReferenceCell(c)),
					1e-8));
			}
		}
	}
	
	@Test
	public void testMappingRotated3D()
	{
		final DistortedCell cell = getRotated3D();
		final CoordinateDenseMatrix rotation = CoordinateDenseMatrix.fromValues(3, 3, 1, 1, 0, -1, 1, 0, 0, 0,
		                                                                        1);
		final CoordinateMatrix rotationInv = rotation.inverse();
		for (int i = 0; i <= 10; i++)
		{
			for (int j = 0; j <= 10; j++)
			{
				for (int k = 0; k <= 10; k++)
				{
					final CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1,
					                                                       k * 0.1);
					assertEquals(cell.transformationGradientFromReferenceCell(c), rotation);
					assertEquals(cell.transformationGradientToReferenceCell(c), rotationInv);
					assertTrue(c.almostEqual(cell.transformPreciseToReferenceCell(
						cell.transformFromReferenceCell(c)), 1e-8));
				}
			}
		}
	}
	
	@Test
	public void testMappingDistorted3D()
	{
		final DistortedCell cell = getDistortedCell3D();
		for (int i = 0; i <= 10; i++)
		{
			for (int j = 0; j <= 10; j++)
			{
				for (int k = 0; k <= 10; k++)
				{
					final CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1,
					                                                       k * 0.1);
					assertTrue(c.almostEqual(cell.transformPreciseToReferenceCell(
						cell.transformFromReferenceCell(c)), 1e-8));
				}
			}
		}
	}
	
	@Test
	public void testAnotherMappingDistorted2D()
	{
		final DistortedCell cell = getAnotherDistortedCell2D();
		for (int i = 0; i <= 10; i++)
		{
			for (int j = 0; j <= 10; j++)
			{
				final CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1);
				assertTrue(c.almostEqual(
					cell.transformPreciseToReferenceCell(cell.transformFromReferenceCell(c)),
					1e-8));
			}
		}
	}
	
	@Test
	public void testMappingDistorted2D()
	{
		final DistortedCell cell = getDistortedCell2D();
		for (int i = 0; i <= 10; i++)
		{
			for (int j = 0; j <= 10; j++)
			{
				final CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1);
				assertTrue(c.almostEqual(
					cell.transformPreciseToReferenceCell(cell.transformFromReferenceCell(c)),
					1e-8));
			}
		}
	}
}
