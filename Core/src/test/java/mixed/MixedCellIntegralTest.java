package mixed;

import basic.ScalarFunction;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.junit.Test;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPCellIntegral;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.TPCell;

import static org.junit.Assert.assertTrue;

public class MixedCellIntegralTest
{
	@Test
	public void testPressureIntegral()
	{
		final CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(3, 0),
		                                          CoordinateVector.fromValues(5, 4), new IntCoordinates(1, 1));
		final TPCell c = g.cells.get(0);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					final ContinuousTPShapeFunction f1 = new ContinuousTPShapeFunction(c, k, i);
					final ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c, k, j);
					final QkQkFunction mf1 = new QkQkFunction(f1);
					final QkQkFunction mf2 = new QkQkFunction(f2);
					final TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(
						ScalarFunction.constantFunction(1),
						TPCellIntegral.GRAD_GRAD);
					final MixedCellIntegral<TPCell, ContinuousTPShapeFunction,
						ContinuousTPVectorFunction, QkQkFunction> gg2 =
						MixedCellIntegral.fromPressureIntegral(gg);
					assertTrue(Math.abs(
						gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						                                                              mf1,
						                                                              mf2)) < 1e-12);
				}
	}
	
	@Test
	public void testVelocityIntegral()
	{
		final CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(3, 0),
		                                          CoordinateVector.fromValues(5, 4), new IntCoordinates(1, 1));
		final TPCell c = g.cells.get(0);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					final ContinuousTPVectorFunction f1 = new ContinuousTPVectorFunction(c, k, i);
					final ContinuousTPVectorFunction f2 = new ContinuousTPVectorFunction(c, k, j);
					final QkQkFunction mf1 = new QkQkFunction(f1);
					final QkQkFunction mf2 = new QkQkFunction(f2);
					final TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
						new TPVectorCellIntegral<>(ScalarFunction.constantFunction(1),
						                           TPCellIntegral.GRAD_GRAD);
					final MixedCellIntegral<TPCell, ContinuousTPShapeFunction,
						ContinuousTPVectorFunction, QkQkFunction> gg2 =
						MixedCellIntegral.fromVelocityIntegral(gg);
					assertTrue(Math.abs(
						gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						                                                              mf1,
						                                                              mf2)) < 1e-12);
				}
	}
	
	@Test
	public void testMixedIntegral()
	{
		final CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(3, 0),
		                                          CoordinateVector.fromValues(5, 4), new IntCoordinates(1, 1));
		final TPCell c = g.cells.get(0);
		
		ContinuousTPVectorFunction f1 = null;
		final MixedTPCellIntegral<ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vg =
			new MixedTPCellIntegral<>(MixedTPCellIntegral.VALUE_GRAD);
		final MixedTPCellIntegral<ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> dv =
			new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		final double[] valuesVG = new double[]{-0.666, -0.666, -0.333, -0.333, -0.333, -0.1666, -0.333, -0.1666};
		final double[] valuesDV = new double[]{-0.666, 0.666, -0.333, 0.333, -0.333, -0.1666, 0.333, 0.1666};
		for (int i = 0; i < 8; i++)
		{
			f1 = new ContinuousTPVectorFunction(c, 1, i);
			final ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c, 1, 0);
			final QkQkFunction mf1 = new QkQkFunction(f1);
			final QkQkFunction mf2 = new QkQkFunction(f2);
			assertTrue(Math.abs(vg.evaluateCellIntegral(c, mf1, mf2) - valuesVG[i]) < 1e-2);
			assertTrue(Math.abs(dv.evaluateCellIntegral(c, mf1, mf2) - valuesDV[i]) < 1e-2);
		}
	}
}
