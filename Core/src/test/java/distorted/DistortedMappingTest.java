package distorted;

import distorted.geometry.DistortedCell;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class DistortedMappingTest
{
	@Test
	public void testMappingNonDistorted3D()
	{
		DistortedCell cell = getNonDistortedCell3D();
		for(int i = 0; i <= 10; i++)
		{
			for(int j = 0; j <= 10; j++)
			{
				for(int k = 0; k <= 10; k++)
				{
					CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1, k*0.1);
					assertTrue(c.almostEqual(cell.transformToReferenceCell(c), 1e-8));
					assertTrue(c.almostEqual(cell.transformFromReferenceCell(c), 1e-8));
					assertEquals(cell.transformationGradientFromReferenceCell(c),
						CoordinateMatrix.fromValues(3, 3,
							1, 0, 0,
							0, 1, 0,
							0, 0, 1));
					assertEquals(cell.transformationGradientToReferenceCell(c),CoordinateMatrix.fromValues(3, 3,
						1, 0, 0,
						0, 1, 0,
						0, 0, 1));
					assertTrue(c.almostEqual(cell.transformToReferenceCell(cell.transformFromReferenceCell(c)), 1e-8));
				}
			}
		}
	}
	@Test
	public void testMappingNonDistorted2D()
	{
		DistortedCell cell = getNonDistortedCell2D();
		for(int i = 0; i <=10; i++)
		{
			for(int j = 0; j <=10; j++)
			{
				CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1);
				assertTrue(c.almostEqual(cell.transformToReferenceCell(c),1e-8));
				assertTrue(c.almostEqual(cell.transformFromReferenceCell(c),1e-8));
				assertEquals(cell.transformationGradientFromReferenceCell(c),
					CoordinateMatrix.fromValues(2,2,1,0,0,1));
				assertEquals(cell.transformationGradientToReferenceCell(c),
					CoordinateMatrix.fromValues(2,2,1,0,0,1));
				assertTrue(c.almostEqual(cell.transformToReferenceCell(cell.transformFromReferenceCell(c)),1e-8));
			}
		}
	}
	@Test
	public void testMappingScaled2D()
	{
		DistortedCell cell = getScaled2D();
		CoordinateMatrix scale = CoordinateMatrix.fromValues(2,2,4,0,0,2);
		CoordinateMatrix scaleInv = CoordinateMatrix.fromValues(2,2,0.25,0,0,0.5);
		for(int i = 0; i <=10; i++)
		{
			for(int j = 0; j <=10; j++)
			{
				CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1);
				assertEquals(cell.transformationGradientFromReferenceCell(c),
					scale);
				assertEquals(cell.transformationGradientToReferenceCell(c),
					scaleInv);
				assertTrue(c.almostEqual(cell.transformToReferenceCell(scale.mvMul(c)),1e-8));
				assertTrue(c.almostEqual(cell.transformFromReferenceCell(scaleInv.mvMul(c)),1e-8));
				assertTrue(c.almostEqual(cell.transformToReferenceCell(cell.transformFromReferenceCell(c)),1e-8));
			}
		}
	}
	@Test
	public void testMappingScaled3D()
	{
		DistortedCell cell = getScaled3D();
		CoordinateMatrix scale = CoordinateMatrix.fromValues(3,3,4,0,0,0,2,0,0,0,8);
		CoordinateMatrix scaleInv = CoordinateMatrix.fromValues(3,3,0.25,0,0,0,0.5,0,0,0,0.125);
		for(int i = 0; i <= 10; i++)
		{
			for(int j = 0; j <= 10; j++)
			{
				for(int k = 0; k <= 10; k++)
				{
					CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1, k*0.1);
					assertTrue(c.almostEqual(cell.transformToReferenceCell(scale.mvMul(c)),1e-8));
					assertTrue(c.almostEqual(cell.transformFromReferenceCell(scaleInv.mvMul(c)),1e-8));
					assertEquals(cell.transformationGradientFromReferenceCell(c),scale);
					assertEquals(cell.transformationGradientToReferenceCell(c),scaleInv);
					assertTrue(c.almostEqual(cell.transformToReferenceCell(cell.transformFromReferenceCell(c)), 1e-8));
				}
			}
		}
	}
	@Test
	public void testMappingRotated2D()
	{
		DistortedCell cell = getRotated2D();
		CoordinateMatrix rotation = CoordinateMatrix.fromValues(2,2,1,1,-1,1);
		CoordinateMatrix rotationInv = rotation.inverse();
		for(int i = 0; i <=10; i++)
		{
			for(int j = 0; j <=10; j++)
			{
				CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1);
				assertEquals(cell.transformationGradientFromReferenceCell(c),
					rotation);
				assertEquals(cell.transformationGradientToReferenceCell(c),
					rotationInv);
				assertTrue(c.almostEqual(cell.transformToReferenceCell(cell.transformFromReferenceCell(c)),1e-8));
			}
		}
	}
	@Test
	public void testMappingRotated3D()
	{
		DistortedCell cell = getRotated3D();
		CoordinateMatrix rotation = CoordinateMatrix.fromValues(3,3,1,1,0,-1,1,0,0,0,1);
		CoordinateMatrix rotationInv = rotation.inverse();
		for(int i = 0; i <= 10; i++)
		{
			for(int j = 0; j <= 10; j++)
			{
				for(int k = 0; k <= 10; k++)
				{
					CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1, k*0.1);
					assertEquals(cell.transformationGradientFromReferenceCell(c), rotation);
					assertEquals(cell.transformationGradientToReferenceCell(c), rotationInv);
					assertTrue(c.almostEqual(cell.transformToReferenceCell(cell.transformFromReferenceCell(c)), 1e-8));
				}
			}
		}
	}
	@Test
	public void testMappingDistorted3D()
	{
		DistortedCell cell = getDistortedCell3D();
		for(int i = 0; i <= 10; i++)
		{
			for(int j = 0; j <= 10; j++)
			{
				for(int k = 0; k <= 10; k++)
				{
					CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1, k*0.1);
					assertTrue(c.almostEqual(cell.transformToReferenceCell(cell.transformFromReferenceCell(c)), 1e-8));
				}
			}
		}
	}
	@Test
	public void testMappingDistorted2D()
	{
		DistortedCell cell = getDistortedCell2D();
		for(int i = 0; i <=10; i++)
		{
			for(int j = 0; j <=10; j++)
			{
				CoordinateVector c = CoordinateVector.fromValues(i * 0.1, j * 0.1);
				assertTrue(c.almostEqual(cell.transformToReferenceCell(cell.transformFromReferenceCell(c)),1e-8));
			}
		}
	}
	
	private DistortedCell getNonDistortedCell2D()
	{
		CoordinateVector [] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0,0);
		vertices[1] = CoordinateVector.fromValues(1,0);
		vertices[2] = CoordinateVector.fromValues(1,1);
		vertices[3] = CoordinateVector.fromValues(0,1);
		return new DistortedCell(vertices);
	}
	private DistortedCell getNonDistortedCell3D()
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
		return new DistortedCell(vertices);
	}
	private DistortedCell getScaled2D()
	{
		CoordinateVector [] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0,0);
		vertices[1] = CoordinateVector.fromValues(4,0);
		vertices[2] = CoordinateVector.fromValues(4,2);
		vertices[3] = CoordinateVector.fromValues(0,2);
		return new DistortedCell(vertices);
	}
	private DistortedCell getScaled3D()
	{
		CoordinateVector [] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0,0,0);
		vertices[1] = CoordinateVector.fromValues(4,0,0);
		vertices[2] = CoordinateVector.fromValues(4,2,0);
		vertices[3] = CoordinateVector.fromValues(0,2,0);
		vertices[4] = CoordinateVector.fromValues(0,0,8);
		vertices[5] = CoordinateVector.fromValues(4,0,8);
		vertices[6] = CoordinateVector.fromValues(4,2,8);
		vertices[7] = CoordinateVector.fromValues(0,2,8);
		return new DistortedCell(vertices);
	}
	private DistortedCell getRotated2D()
	{
		CoordinateVector [] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0,0);
		vertices[1] = CoordinateVector.fromValues(1,-1);
		vertices[2] = CoordinateVector.fromValues(2,0);
		vertices[3] = CoordinateVector.fromValues(1,1);
		return new DistortedCell(vertices);
	}
	private DistortedCell getRotated3D()
	{
		CoordinateVector [] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(0,0,0);
		vertices[1] = CoordinateVector.fromValues(1,-1,0);
		vertices[2] = CoordinateVector.fromValues(2,0,0);
		vertices[3] = CoordinateVector.fromValues(1,1,0);
		vertices[4] = CoordinateVector.fromValues(0,0,1);
		vertices[5] = CoordinateVector.fromValues(1,-1,1);
		vertices[6] = CoordinateVector.fromValues(2,0,1);
		vertices[7] = CoordinateVector.fromValues(1,1,1);
		return new DistortedCell(vertices);
	}
	private DistortedCell getDistortedCell2D()
	{
		CoordinateVector [] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0,0);
		vertices[1] = CoordinateVector.fromValues(2,0);
		vertices[2] = CoordinateVector.fromValues(0.8,23.1);
		vertices[3] = CoordinateVector.fromValues(0,1);
		return new DistortedCell(vertices);
	}
	private DistortedCell getDistortedCell3D()
	{
		CoordinateVector [] vertices = new CoordinateVector[8];
		vertices[0] = CoordinateVector.fromValues(5,5,5);
		vertices[1] = CoordinateVector.fromValues(10,0,5);
		vertices[2] = CoordinateVector.fromValues(15,5,5);
		vertices[3] = CoordinateVector.fromValues(10,10,5);
		vertices[4] = CoordinateVector.fromValues(4,6,11);
		vertices[5] = CoordinateVector.fromValues(9,1,11);
		vertices[6] = CoordinateVector.fromValues(14,6,11);
		vertices[7] = CoordinateVector.fromValues(9,11,11);
		return new DistortedCell(vertices);
	}
}
