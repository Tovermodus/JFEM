package tensorproduct;

import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import org.junit.Test;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import static org.junit.Assert.*;

public class TPFESpaceTest
{
	
	@Test
	public void testSquareSpace()
	{
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 1;
		final TPFESpace grid = new TPFESpace(start, end,
		                                     Ints.asList(3, 3));
		grid.assembleCells();
		;
		grid.assembleFunctions(polynomialDegree);
		assertEquals(9, grid.getCells().size());
		assertEquals(24, grid.getFaces().size());
		assertEquals(9 * (polynomialDegree + 1) * (polynomialDegree + 1), grid.shapeFunctions.size());
		for (final TPCell cell : grid.getCells())
		{
			assertEquals(grid.getShapeFunctionsWithSupportOnCell(cell).size(),
			             (polynomialDegree + 1) * (polynomialDegree + 1));
			for (final TPShapeFunction shapeFunction : grid.getShapeFunctionsWithSupportOnCell(cell))
			{
				assertTrue(cell.isInCell(shapeFunction.nodeFunctional.getPoint()));
			}
		}
		for (final TPCell cell : grid.getCells())
		{
			for (final TPFace face : grid.getFaces())
			{
				if (cell.isInCell(face.center()))
				{
					assertTrue(cell.getFaces().contains(face));
					assertTrue(face.getCells().contains(cell));
				} else
				{
					assertFalse(cell.getFaces().contains(face));
					assertFalse(face.getCells().contains(cell));
				}
			}
			for (final TPShapeFunction shapeFunction : grid.getShapeFunctions().values())
			{
				if (shapeFunction.value(cell
					                        .center()
					                        .add(CoordinateVector.fromValues(Math.random() * 1e-4,
					                                                         Math.random() * 1e-4))) != 0)
					assertTrue(
						grid.getShapeFunctionsWithSupportOnCell(cell).contains(shapeFunction));
			}
		}
		for (final TPFace face : grid.getFaces())
		{
			
			for (final TPShapeFunction shapeFunction : grid.getShapeFunctions().values())
			{
				if (shapeFunction.value(
					face.center().add(CoordinateVector.fromValues(Math.random() * 1e-8,
					                                              Math.random() * 1e-8))) != 0)
					assertTrue(grid.getShapeFunctionsWithSupportOnFace(
						face
					).contains(shapeFunction));
			}
			if (face.isBoundaryFace())
				assertEquals(grid.getShapeFunctionsWithSupportOnFace(
					face
				).size(), (polynomialDegree + 1) * (polynomialDegree + 1));
			else
				assertEquals(grid.getShapeFunctionsWithSupportOnFace(
					face
				).size(), (polynomialDegree + 1) * (polynomialDegree + 1) * 2);
		}
	}
	
	@Test
	public void testSquareSpaceDeg3()
	{
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 3;
		final TPFESpace grid = new TPFESpace(start, end,
		                                     Ints.asList(3, 3));
		grid.assembleCells();
		;
		grid.assembleFunctions(polynomialDegree);
		assertEquals(9, grid.getCells().size());
		assertEquals(24, grid.getFaces().size());
		assertEquals(9 * (polynomialDegree + 1) * (polynomialDegree + 1), grid.shapeFunctions.size());
		for (final TPCell cell : grid.getCells())
		{
			assertEquals(grid.getShapeFunctionsWithSupportOnCell(cell).size(),
			             (polynomialDegree + 1) * (polynomialDegree + 1));
			for (final TPShapeFunction shapeFunction : grid.getShapeFunctionsWithSupportOnCell(cell))
			{
				assertTrue(cell.isInCell(shapeFunction.nodeFunctional.getPoint()));
			}
		}
		for (final TPCell cell : grid.getCells())
		{
			for (final TPFace face : grid.getFaces())
			{
				if (cell.isInCell(face.center()))
				{
					assertTrue(cell.getFaces().contains(face));
					assertTrue(face.getCells().contains(cell));
				} else
				{
					assertFalse(cell.getFaces().contains(face));
					assertFalse(face.getCells().contains(cell));
				}
			}
			for (final TPShapeFunction shapeFunction : grid.getShapeFunctions().values())
			{
				if (shapeFunction.value(cell
					                        .center()
					                        .add(CoordinateVector.fromValues(Math.random() * 1e-4,
					                                                         Math.random() * 1e-4))) != 0)
					assertTrue(
						grid.getShapeFunctionsWithSupportOnCell(cell).contains(shapeFunction));
			}
		}
		for (final TPFace face : grid.getFaces())
		{
			
			for (final TPShapeFunction shapeFunction : grid.getShapeFunctions().values())
			{
				if (shapeFunction.value(
					face.center().add(CoordinateVector.fromValues(Math.random() * 1e-8,
					                                              Math.random() * 1e-8))) != 0)
					assertTrue(grid.getShapeFunctionsWithSupportOnFace(
						face
					).contains(shapeFunction));
			}
			if (face.isBoundaryFace())
				assertEquals(grid.getShapeFunctionsWithSupportOnFace(
					face
				).size(), (polynomialDegree + 1) * (polynomialDegree + 1));
			else
				assertEquals(grid.getShapeFunctionsWithSupportOnFace(
					face
				).size(), (polynomialDegree + 1) * (polynomialDegree + 1) * 2);
		}
	}
}
