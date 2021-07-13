package mixed;

import basic.PlotFrame;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;

import org.junit.jupiter.api.Test;

import java.util.Map;
import java.util.TreeMap;

import static org.junit.jupiter.api.Assertions.*;

public class NedelecSpaceTest
{
	@Test
	public void testOneCell()
	{
		CoordinateVector start = CoordinateVector.fromValues(0, 0,0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1,1);
		int polynomialDegree = 2;
		NedelecSpace grid = new NedelecSpace(start, end,
			Ints.asList(1,1,1), polynomialDegree);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		assertEquals(3*polynomialDegree*(polynomialDegree+1)*(polynomialDegree+1)+Math.pow(polynomialDegree
			,3),
			grid.shapeFunctions.size());
		assertEquals(3*polynomialDegree*(polynomialDegree+1)*(polynomialDegree+1),grid.nodeFuncionals.size());
		assertEquals(1,grid.getCells().size());
		assertEquals(6,grid.getFaces().size());
		assertEquals(12,grid.edges.size());
		for(MixedShapeFunction<?,?,?,?,?> shapeFunction: grid.shapeFunctions)
		{
			if(shapeFunction.isVelocity())
			{
				for (NedelecNodeFuncional n : grid.nodeFuncionals)
				{
					if(grid.functionalShapeFunctionMap.get(n).equals(shapeFunction.getVelocityShapeFunction()))
						assertEquals(1.0,
							n.evaluate
								(shapeFunction.getVelocityShapeFunction()), 1e-12);
					else
						assertEquals(0.0,
							n.evaluate
								(shapeFunction.getVelocityShapeFunction()), 1e-12);
				}
			}
			assertEquals(1.0, shapeFunction.getNodeFunctional().evaluate(shapeFunction),1e-12);
		}
		
	}
	@Test
	public void testTwoCells()
	{
		CoordinateVector start = CoordinateVector.fromValues(0, 0,0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1,1);
		int polynomialDegree = 2;
		NedelecSpace grid = new NedelecSpace(start, end,
			Ints.asList(2,2,2), polynomialDegree);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		//assertEquals(2,grid.getCells().size());
		//assertEquals(11,grid.getFaces().size());
		//assertEquals(20,grid.edges.size());
		for(MixedShapeFunction<?,?,?,?,?> shapeFunction: grid.shapeFunctions)
		{
			if(shapeFunction.isVelocity())
			{
				for (NedelecNodeFuncional n : grid.nodeFuncionals)
				{
					if(grid.functionalShapeFunctionMap.get(n).equals(shapeFunction.getVelocityShapeFunction()))
						assertEquals(1.0,
							n.evaluate
								(shapeFunction.getVelocityShapeFunction()), 1e-12);
					else
						assertEquals(0.0,
							n.evaluate
								(shapeFunction.getVelocityShapeFunction()), 1e-12);
				}
			}
			assertEquals(1.0,shapeFunction.getNodeFunctional().evaluate(shapeFunction),1e-12);
		}
		
	}
}
