package distorted;

import distorted.geometry.DistortedCell;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.DenseVector;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class DistortedCellTest
{
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
	private DistortedCell getScaled2D()
	{
		CoordinateVector [] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(0,0);
		vertices[1] = CoordinateVector.fromValues(4,0);
		vertices[2] = CoordinateVector.fromValues(4,2);
		vertices[3] = CoordinateVector.fromValues(0,2);
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
}
