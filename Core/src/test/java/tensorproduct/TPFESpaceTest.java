package tensorproduct;

import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import org.junit.jupiter.api.Test;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import static org.junit.jupiter.api.Assertions.*;

public class TPFESpaceTest
{
	
	@Test
	public void testSquareSpace()
	{
		CoordinateVector start = CoordinateVector.fromValues(0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 1;
		TPFESpace grid = new TPFESpace(start, end,
			Ints.asList(3, 3));
		grid.assembleCells();;
		grid.assembleFunctions(polynomialDegree);
		assertEquals(9, grid.getCells().size());
		assertEquals(24, grid.getFaces().size());
		assertEquals(9*(polynomialDegree+1)*(polynomialDegree+1), grid.shapeFunctions.size());
		for(TPCell cell: grid.getCells())
		{
			assertEquals(grid.getShapeFunctionsWithSupportOnCell(cell).size(),
				(polynomialDegree + 1) * (polynomialDegree + 1));
			for(TPShapeFunction shapeFunction: grid.getShapeFunctionsWithSupportOnCell(cell))
			{
				assertTrue(cell.isInCell(shapeFunction.nodeFunctional.getPoint()));
			}
		}
		for(TPCell cell: grid.getCells())
		{
			for (TPFace face: grid.getFaces())
			{
				if(cell.isInCell(face.center()))
				{
					assertTrue(cell.getFaces().contains(face));
					assertTrue(face.getCells().contains(cell));
				}
				else
				{
					assertFalse(cell.getFaces().contains(face));
					assertFalse(face.getCells().contains(cell));
				}
			}
			for(TPShapeFunction shapeFunction: grid.getShapeFunctions().values())
			{
				if(shapeFunction.value(cell.center().add(CoordinateVector.fromValues(Math.random()*1e-4, Math.random()*1e-4))) != 0)
					assertTrue(grid.getShapeFunctionsWithSupportOnCell(cell).contains(shapeFunction));
			}
		}
		for(TPFace face:grid.getFaces())
		{
			
			for(TPShapeFunction shapeFunction: grid.getShapeFunctions().values())
			{
				if(shapeFunction.value(face.center().add(CoordinateVector.fromValues(Math.random()*1e-8,
					Math.random()*1e-8))) != 0)
					assertTrue(grid.getShapeFunctionsWithSupportOnFace(
						face
					).contains(shapeFunction));
			}
			if(face.isBoundaryFace())
				assertEquals(grid.getShapeFunctionsWithSupportOnFace(
					face
				).size(), (polynomialDegree+1)*(polynomialDegree+1));
			else
				assertEquals(grid.getShapeFunctionsWithSupportOnFace(
					face
				).size(), (polynomialDegree+1)*(polynomialDegree+1)*2);
			
		}
		
		
	}
	@Test
	public void testSquareSpaceDeg3()
	{
		CoordinateVector start = CoordinateVector.fromValues(0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 3;
		TPFESpace grid = new TPFESpace(start, end,
			Ints.asList(3, 3));
		grid.assembleCells();;
		grid.assembleFunctions(polynomialDegree);
		assertEquals(9, grid.getCells().size());
		assertEquals(24, grid.getFaces().size());
		assertEquals(9*(polynomialDegree+1)*(polynomialDegree+1), grid.shapeFunctions.size());
		for(TPCell cell: grid.getCells())
		{
			assertEquals(grid.getShapeFunctionsWithSupportOnCell(cell).size(),
				(polynomialDegree + 1) * (polynomialDegree + 1));
			for(TPShapeFunction shapeFunction: grid.getShapeFunctionsWithSupportOnCell(cell))
			{
				assertTrue(cell.isInCell(shapeFunction.nodeFunctional.getPoint()));
			}
		}
		for(TPCell cell: grid.getCells())
		{
			for (TPFace face: grid.getFaces())
			{
				if(cell.isInCell(face.center()))
				{
					assertTrue(cell.getFaces().contains(face));
					assertTrue(face.getCells().contains(cell));
				}
				else
				{
					assertFalse(cell.getFaces().contains(face));
					assertFalse(face.getCells().contains(cell));
				}
			}
			for(TPShapeFunction shapeFunction: grid.getShapeFunctions().values())
			{
				if(shapeFunction.value(cell.center().add(CoordinateVector.fromValues(Math.random()*1e-4, Math.random()*1e-4))) != 0)
					assertTrue(grid.getShapeFunctionsWithSupportOnCell(cell).contains(shapeFunction));
			}
		}
		for(TPFace face:grid.getFaces())
		{
			
			for(TPShapeFunction shapeFunction: grid.getShapeFunctions().values())
			{
				if(shapeFunction.value(face.center().add(CoordinateVector.fromValues(Math.random()*1e-8,
					Math.random()*1e-8))) != 0)
					assertTrue(grid.getShapeFunctionsWithSupportOnFace(
						face
					).contains(shapeFunction));
			}
			if(face.isBoundaryFace())
				assertEquals(grid.getShapeFunctionsWithSupportOnFace(
					face
				).size(), (polynomialDegree+1)*(polynomialDegree+1));
			else
				assertEquals(grid.getShapeFunctionsWithSupportOnFace(
					face
				).size(), (polynomialDegree+1)*(polynomialDegree+1)*2);
			
		}
		
		
	}
	
	
	
}