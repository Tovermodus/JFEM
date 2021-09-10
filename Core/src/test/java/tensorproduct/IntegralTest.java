package tensorproduct;

import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.junit.Test;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.TPCell;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class IntegralTest
{
	
	@Test
	public void testTPCellIntegralDeg1G2()
	{
		final int polynomialDegree = 1;
		final double[] testValuesgg = new double[]{
			0.6833, -0.0583, -0.2833, -0.3417, -0.0583, 0.6833, -0.3417, -0.2833, -0.2833, -0.3417,
			0.6833, -0.0583, -0.3417, -0.2833, -0.0583, 0.6833};
		final double[] testValuesggw = new double[]{
			7.7812, -0.4688, -3.1612, -4.1513, -0.4688, 8.3437, -4.1513, -3.7237, -3.1612, -4.1513,
			8.2612, -0.9488, -4.1513, -3.7237, -0.9488, 8.8238};
		final double[] testValuesvv = new double[]{
			0.0889, 0.0444, 0.0444, 0.0222, 0.0444, 0.0889, 0.0222, 0.0444,
			0.0444, 0.0222, 0.0889, 0.0444, 0.0222, 0.0444, 0.0444, 0.0889};
		final double[] testValuesvvw = new double[]{
			0.9444, 0.5000, 0.5100, 0.2700, 0.5000, 1.0556, 0.2700, 0.5700, 0.5100, 0.2700, 1.0956,
			0.5800, 0.2700, 0.5700, 0.5800, 1.2244};
		final double[] testValuesvg = new double[]{
			-0.3625, 0.1486, -0.0222, 0.2361, -0.3181, -0.1097, 0.0001, 0.4277, -0.2472, 0.0111, -0.1428,
			0.3789, -0.2249, -0.2473, -0.0984, 0.5706};
		final double[] testValuesgv = new double[]{
			-0.3625, -0.3181, -0.2472, -0.2249, 0.1486, -0.1097, 0.0111, -0.2473, -0.0222, 0.0001,
			-0.1428, -0.0984, 0.2361, 0.4277, 0.3789, 0.5706};
		
		final CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0, 2.3),
		                                          CoordinateVector.fromValues(1, 3.1), new IntCoordinates(1, 1));
		final TPCell cell = g.cells.get(0);
		final List<TPShapeFunction> shapeFunctions = new ArrayList<>();
		for (int i = 0; i < Math.pow(polynomialDegree + 1, 2); i++)
		{
			shapeFunctions.add(new TPShapeFunction(cell, polynomialDegree, i));
		}
		final ScalarFunction weight = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return pos.at(1) * (pos.at(0) + 4);
			}
		};
		final VectorFunction vweight = new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(2 - pos.x() + pos.y() * 0.1, pos.x() * pos.y());
			}
		};
		final TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD,
		                                                                QuadratureRule1D.Gauss2);
		final TPCellIntegral<TPShapeFunction> ggw = new TPCellIntegral<>(weight, TPCellIntegral.GRAD_GRAD,
		                                                                 QuadratureRule1D.Gauss2);
		final TPCellIntegral<TPShapeFunction> vg = new TPCellIntegral<>(vweight, TPCellIntegral.VALUE_GRAD,
		                                                                QuadratureRule1D.Gauss2);
		final TPCellIntegral<TPShapeFunction> gv = new TPCellIntegral<>(vweight, TPCellIntegral.GRAD_VALUE,
		                                                                QuadratureRule1D.Gauss2);
		final TPCellIntegral<TPShapeFunction> vv = new TPCellIntegral<>(TPCellIntegral.VALUE_VALUE,
		                                                                QuadratureRule1D.Gauss2);
		final TPCellIntegral<TPShapeFunction> vvw = new TPCellIntegral<>(weight, TPCellIntegral.VALUE_VALUE,
		                                                                 QuadratureRule1D.Gauss2);
		int k = 0;
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				
				assertTrue(Math.abs(gg.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                            shapeFunctions.get(
					                                            j)) - testValuesgg[k]) <= 1e-2);
				assertTrue(Math.abs(ggw.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                             shapeFunctions.get(
					                                             j)) - testValuesggw[k]) <= 1e-2);
				assertTrue(Math.abs(vg.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                            shapeFunctions.get(
					                                            j)) - testValuesvg[k]) <= 1e-2);
				assertTrue(Math.abs(gv.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                            shapeFunctions.get(
					                                            j)) - testValuesgv[k]) <= 1e-2);
				assertTrue(Math.abs(vv.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                            shapeFunctions.get(
					                                            j)) - testValuesvv[k]) <= 1e-2);
				assertTrue(Math.abs(vvw.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                             shapeFunctions.get(
					                                             j)) - testValuesvvw[k++]) <= 1e-2);
			}
		}
	}
	
	@Test
	public void testTPCellIntegralDeg1G3()
	{
		final int polynomialDegree = 1;
		final double[] testValuesgg = new double[]{
			0.6833, -0.0583, -0.2833, -0.3417, -0.0583, 0.6833, -0.3417, -0.2833, -0.2833, -0.3417,
			0.6833, -0.0583, -0.3417, -0.2833, -0.0583, 0.6833};
		final double[] testValuesggw = new double[]{
			7.7812, -0.4688, -3.1612, -4.1513, -0.4688, 8.3437, -4.1513, -3.7237, -3.1612, -4.1513,
			8.2612, -0.9488, -4.1513, -3.7237, -0.9488, 8.8238};
		final double[] testValuesvv = new double[]{
			0.0889, 0.0444, 0.0444, 0.0222, 0.0444, 0.0889, 0.0222, 0.0444,
			0.0444, 0.0222, 0.0889, 0.0444, 0.0222, 0.0444, 0.0444, 0.0889};
		final double[] testValuesvvw = new double[]{
			0.9444, 0.5000, 0.5100, 0.2700, 0.5000, 1.0556, 0.2700, 0.5700, 0.5100, 0.2700, 1.0956,
			0.5800, 0.2700, 0.5700, 0.5800, 1.2244};
		final double[] testValuesvg = new double[]{
			-0.3625, 0.1486, -0.0222, 0.2361, -0.3181, -0.1097, 0.0001, 0.4277, -0.2472, 0.0111, -0.1428,
			0.3789, -0.2249, -0.2473, -0.0984, 0.5706};
		final double[] testValuesgv = new double[]{
			-0.3625, -0.3181, -0.2472, -0.2249, 0.1486, -0.1097, 0.0111, -0.2473, -0.0222, 0.0001,
			-0.1428, -0.0984, 0.2361, 0.4277, 0.3789, 0.5706};
		
		final CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(0, 2.3),
		                                          CoordinateVector.fromValues(1, 3.1), new IntCoordinates(1, 1));
		final TPCell cell = g.cells.get(0);
		final List<TPShapeFunction> shapeFunctions = new ArrayList<>();
		for (int i = 0; i < Math.pow(polynomialDegree + 1, 2); i++)
		{
			shapeFunctions.add(new TPShapeFunction(cell, polynomialDegree, i));
		}
		final ScalarFunction weight = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return pos.at(1) * (pos.at(0) + 4);
			}
		};
		final VectorFunction vweight = new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(2 - pos.x() + pos.y() * 0.1, pos.x() * pos.y());
			}
		};
		final TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD,
		                                                                QuadratureRule1D.Gauss3);
		final TPCellIntegral<TPShapeFunction> ggw = new TPCellIntegral<>(weight, TPCellIntegral.GRAD_GRAD,
		                                                                 QuadratureRule1D.Gauss3);
		final TPCellIntegral<TPShapeFunction> vg = new TPCellIntegral<>(vweight, TPCellIntegral.VALUE_GRAD,
		                                                                QuadratureRule1D.Gauss3);
		final TPCellIntegral<TPShapeFunction> gv = new TPCellIntegral<>(vweight, TPCellIntegral.GRAD_VALUE,
		                                                                QuadratureRule1D.Gauss3);
		final TPCellIntegral<TPShapeFunction> vv = new TPCellIntegral<>(TPCellIntegral.VALUE_VALUE,
		                                                                QuadratureRule1D.Gauss3);
		final TPCellIntegral<TPShapeFunction> vvw = new TPCellIntegral<>(weight, TPCellIntegral.VALUE_VALUE,
		                                                                 QuadratureRule1D.Gauss3);
		int k = 0;
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				
				assertTrue(Math.abs(gg.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                            shapeFunctions.get(
					                                            j)) - testValuesgg[k]) <= 1e-2);
				assertTrue(Math.abs(ggw.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                             shapeFunctions.get(
					                                             j)) - testValuesggw[k]) <= 1e-2);
				assertTrue(Math.abs(vg.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                            shapeFunctions.get(
					                                            j)) - testValuesvg[k]) <= 1e-2);
				assertTrue(Math.abs(gv.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                            shapeFunctions.get(
					                                            j)) - testValuesgv[k]) <= 1e-2);
				assertTrue(Math.abs(vv.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                            shapeFunctions.get(
					                                            j)) - testValuesvv[k]) <= 1e-2);
				assertTrue(Math.abs(vvw.evaluateCellIntegral(cell, shapeFunctions.get(i),
				                                             shapeFunctions.get(
					                                             j)) - testValuesvvw[k++]) <= 1e-2);
			}
		}
	}
}
