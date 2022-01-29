package distorted;

import basic.LagrangeNodeFunctional;
import basic.ScalarFunction;
import distorted.geometry.DistortedCell;
import kotlin.Pair;
import linalg.CoordinateVector;
import linalg.Vector;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class DistortedShapeFunctionTest
{
	private static ScalarFunction getFunctionOnGrid()
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return 1 - (pos.x() - 5) / 2 + (pos.x() - 5) * (pos.y() - 5) / 4 -
					(pos.y() - 5) * (pos.y() - 5) / 4;
			}
			
			@Override
			public CoordinateVector gradient(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(-0.5 + (pos.y() - 5) / 4,
				                                   (pos.x() - 5) / 4 - (pos.y() - 5) / 2);
			}
		};
	}
	
	private static ScalarFunction getFunctionOnReferenceCell()
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return (1 - pos.x()) * (1. - pos.y());
			}
			
			@Override
			public CoordinateVector gradient(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(pos.y() - 1, pos.x() - 1);
			}
		};
	}
	
	private static Pair<DistortedShapeFunction, DistortedCell> createShapeFunction()
	{
		final CoordinateVector[] vertices = new CoordinateVector[4];
		vertices[0] = CoordinateVector.fromValues(5, 5);
		vertices[1] = CoordinateVector.fromValues(7, 5);
		vertices[2] = CoordinateVector.fromValues(9, 7);
		vertices[3] = CoordinateVector.fromValues(7, 7);
		final DistortedCell cell = new DistortedCell(vertices);
		return new Pair<>(new DistortedShapeFunction(cell, 1, 0), cell);
	}
	
	@Test
	public void testValueInCell()
	{
		final DistortedShapeFunction sf = createShapeFunction().getFirst();
		final DistortedCell cell = createShapeFunction().getSecond();
		
		assertTrue(Math.abs(1 - sf.valueInCell(CoordinateVector.fromValues(5, 5), cell)) < 1e-3);
		assertTrue(Math.abs(0.5 - sf.valueInCell(CoordinateVector.fromValues(6, 5), cell)) < 1e-3);
		assertTrue(Math.abs(0.5 - sf.valueInCell(CoordinateVector.fromValues(6, 6), cell)) < 1e-3);
		assertTrue(Math.abs(0.25 - sf.valueInCell(CoordinateVector.fromValues(7, 6), cell)) < 1e-3);
		assertTrue(Math.abs(0 - sf.valueInCell(CoordinateVector.fromValues(7, 5), cell)) < 1e-3);
		assertTrue(Math.abs(0 - sf.valueInCell(CoordinateVector.fromValues(7, 7), cell)) < 1e-3);
		assertTrue(Math.abs(0 - sf.valueInCell(CoordinateVector.fromValues(9, 7), cell)) < 1e-3);
		for (int i = 0; i < 11; i++)
			for (int j = 0; j < 11; j++)
			{
				final CoordinateVector pos = CoordinateVector.fromValues(5 + 0.2 * (j + i),
				                                                         5 + 0.2 * j);
				assertTrue(Math.abs(getFunctionOnGrid().value(pos) -
					                    sf.valueInCell(pos, cell)) < 1e-3);
			}
	}
	
	@Test
	public void testGradientInCell()
	{
		final DistortedShapeFunction sf = createShapeFunction().getFirst();
		final DistortedCell cell = createShapeFunction().getSecond();
		for (int i = 0; i < 11; i++)
			for (int j = 0; j < 11; j++)
			{
				final CoordinateVector pos = CoordinateVector.fromValues(5 + 0.2 * (j + i),
				                                                         5 + 0.2 * j);
				assertTrue(getFunctionOnGrid().gradient(pos)
				                              .sub(
					                              sf.gradientInCell(pos, cell))
				                              .euclidianNorm() < 1e-3);
			}
	}
	
	@Test
	public void testValueInReferenceCell()
	{
		final DistortedShapeFunction sf = createShapeFunction().getFirst();
		final DistortedCell cell = createShapeFunction().getSecond();
		assertTrue(Math.abs(1 - sf.valueOnReferenceCell(CoordinateVector.fromValues(0, 0), cell)) < 1e-13);
		assertTrue(Math.abs(0.5 - sf.valueOnReferenceCell(CoordinateVector.fromValues(0.5, 0), cell)) < 1e-13);
		assertTrue(Math.abs(0.5 - sf.valueOnReferenceCell(CoordinateVector.fromValues(0, 0.5), cell)) < 1e-13);
		assertTrue(
			Math.abs(0.25 - sf.valueOnReferenceCell(CoordinateVector.fromValues(0.5, 0.5), cell)) < 1e-13);
		assertTrue(Math.abs(0 - sf.valueOnReferenceCell(CoordinateVector.fromValues(1, 0), cell)) < 1e-13);
		assertTrue(Math.abs(0 - sf.valueOnReferenceCell(CoordinateVector.fromValues(0, 1), cell)) < 1e-13);
		assertTrue(Math.abs(0 - sf.valueOnReferenceCell(CoordinateVector.fromValues(1, 1), cell)) < 1e-13);
		for (int i = 0; i < 11; i++)
			for (int j = 0; j < 11; j++)
			{
				final CoordinateVector pos = CoordinateVector.fromValues(0.1 * i, 0.1 * j);
				assertTrue(Math.abs(getFunctionOnReferenceCell().value(pos) -
					                    sf.valueOnReferenceCell(pos, cell)) < 1e-13);
			}
	}
	
	@Test
	public void testContinuity()
	{
	}
	
	@Test
	public void testNodeFunctional()
	{
		final DistortedShapeFunction sf = createShapeFunction().getFirst();
		assertEquals(sf.getNodeFunctional(), new LagrangeNodeFunctional(CoordinateVector.fromValues(5, 5)));
	}
	
	@Test
	public void testNodeFunctionalOnReferenceCell()
	{
		final DistortedShapeFunction sf = createShapeFunction().getFirst();
		final DistortedCell cell = createShapeFunction().getSecond();
		assertEquals(sf.nodeFunctionalOnReferenceCell(cell),
		             new LagrangeNodeFunctional(CoordinateVector.fromValues(0, 0)));
	}
	
	@Test
	public void testCompareDegree1()
	{
		for (int i = 0; i < 3; i++)
		{
			final CircleSpace d = new CircleSpace(CoordinateVector.fromValues(17, 19),
			                                      7,
			                                      i);
			d.assembleFunctions(1);
			final List<CoordinateVector> functionalPoints =
				d
					.getShapeFunctionMap()
					.values()
					.stream()
					.map(fun -> fun.getNodeFunctional()
					               .getPoint())
					.collect(
						Collectors.toList());
			final int finalI = i;
			assertTrue(functionalPoints.stream()
			                           .parallel()
			                           .mapToDouble(point -> functionalPoints.stream()
			                                                                 .filter(other -> !other.equals(
				                                                                 point))
			                                                                 .map(point::sub)
			                                                                 .mapToDouble(
				                                                                 Vector::euclidianNorm)
			                                                                 .min()
			                                                                 .orElse(1))
			                           .allMatch(distance -> distance > Math.pow(0.5, finalI)));
		}
	}
	
	@Test
	public void testCompareDegree2()
	{
		for (int i = 0; i < 3; i++)
		{
			final CircleSpace d = new CircleSpace(CoordinateVector.fromValues(17, -19),
			                                      0.07,
			                                      i);
			d.assembleFunctions(2);
			final List<CoordinateVector> functionalPoints =
				d
					.getShapeFunctionMap()
					.values()
					.stream()
					.map(fun -> fun.getNodeFunctional()
					               .getPoint())
					.collect(
						Collectors.toList());
			final int finalI = i;
			assertTrue(functionalPoints.stream()
			                           .parallel()
			                           .mapToDouble(point -> functionalPoints.stream()
			                                                                 .filter(other -> !other.equals(
				                                                                 point))
			                                                                 .map(point::sub)
			                                                                 .mapToDouble(
				                                                                 Vector::euclidianNorm)
			                                                                 .min()
			                                                                 .orElse(1))
			                           .allMatch(distance -> distance > 0.01 * Math.pow(0.5, finalI)));
		}
	}
	
	@Test
	public void testCompareDegree3()
	{
		for (int i = 0; i < 3; i++)
		{
			final CircleSpace d = new CircleSpace(CoordinateVector.fromValues(-17, -19),
			                                      0.7,
			                                      i);
			d.assembleFunctions(3);
			final List<CoordinateVector> functionalPoints =
				d
					.getShapeFunctionMap()
					.values()
					.stream()
					.map(fun -> fun.getNodeFunctional()
					               .getPoint())
					.collect(
						Collectors.toList());
			final int finalI = i;
			assertTrue(functionalPoints.stream()
			                           .parallel()
			                           .mapToDouble(point -> functionalPoints.stream()
			                                                                 .filter(other -> !other.equals(
				                                                                 point))
			                                                                 .map(point::sub)
			                                                                 .mapToDouble(
				                                                                 Vector::euclidianNorm)
			                                                                 .min()
			                                                                 .orElse(1))
			                           .allMatch(distance -> distance > 0.01 * Math.pow(0.5, finalI)));
		}
	}
}
