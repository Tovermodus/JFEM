package distorted;

import basic.PlotWindow;
import distorted.geometry.CircleGrid;
import linalg.CoordinateVector;
import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.TreeSet;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class DistortedShapeFunctionTest
{
	@Test
	public void testCount2DCentral()
	{
		for (int j = 0; j < 4; j++)
		{
			final CircleGrid circle = new CircleGrid(new CoordinateVector(2), 1,
			                                         j);
			System.out.println(circle.faces.size());
			System.out.println(circle.cells.size());
			//System.out.println(circle.faces);
			final HashSet<DistortedShapeFunction> shapeFunctionSet = new HashSet<>();
			circle.cells.stream().parallel().forEach(c ->
			                                         {
				                                         for (int i = 0; i < 4; i++)
					                                         shapeFunctionSet.add(
						                                         new DistortedShapeFunction(c,
						                                                                    1,
						                                                                    i));
			                                         });
			final TreeSet<DistortedShapeFunction> sortedFunctions = new TreeSet<>(shapeFunctionSet);
			assertEquals(sortedFunctions.size(), shapeFunctionSet.size());
		}
	}
	
	@Test
	public void testCount2DCentralScaled()
	{
		for (int j = 0; j < 4; j++)
		{
			final CircleGrid circle = new CircleGrid( new CoordinateVector(2), 1354,
			                                         j);
			System.out.println(circle.faces.size());
			System.out.println(circle.cells.size());
			//System.out.println(circle.faces);
			final HashSet<DistortedShapeFunction> shapeFunctionSet = new HashSet<>();
			circle.cells.stream().parallel().forEach(c ->
			                                         {
				                                         for (int i = 0; i < 4; i++)
					                                         shapeFunctionSet.add(
						                                         new DistortedShapeFunction(c,
						                                                                    1,
						                                                                    i));
			                                         });
			final TreeSet<DistortedShapeFunction> sortedFunctions = new TreeSet<>(shapeFunctionSet);
			assertEquals(sortedFunctions.size(), shapeFunctionSet.size());
		}
	}
	
	@Test
	public void testCount2DMoved()
	{
		for (int j = 0; j < 4; j++)
		{
			final CircleGrid circle = new CircleGrid( CoordinateVector.fromValues(178, -172), 1,
			                                         j);
			System.out.println(circle.faces.size());
			System.out.println(circle.cells.size());
			//System.out.println(circle.faces);
			final HashSet<DistortedShapeFunction> shapeFunctionSet = new HashSet<>();
			circle.cells.stream().parallel().forEach(c ->
			                                         {
				                                         for (int i = 0; i < 4; i++)
					                                         shapeFunctionSet.add(
						                                         new DistortedShapeFunction(c,
						                                                                    1,
						                                                                    i));
			                                         });
			final TreeSet<DistortedShapeFunction> sortedFunctions = new TreeSet<>(shapeFunctionSet);
			assertEquals(sortedFunctions.size(), shapeFunctionSet.size());
		}
	}
	
	@Test
	public void testCount2DMovedScaled()
	{
		for (int j = 0; j < 4; j++)
		{
			final CircleGrid circle = new CircleGrid( CoordinateVector.fromValues(-178, -172), 1312,
			                                         j);
			System.out.println(circle.faces.size());
			System.out.println(circle.cells.size());
			//System.out.println(circle.faces);
			final HashSet<DistortedShapeFunction> shapeFunctionSet = new HashSet<>();
			circle.cells.stream().parallel().forEach(c ->
			                                         {
				                                         for (int i = 0; i < 4; i++)
					                                         shapeFunctionSet.add(
						                                         new DistortedShapeFunction(c,
						                                                                    1,
						                                                                    i));
			                                         });
			final TreeSet<DistortedShapeFunction> sortedFunctions = new TreeSet<>(shapeFunctionSet);
			assertEquals(sortedFunctions.size(), shapeFunctionSet.size());
		}
	}
}
