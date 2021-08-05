package mixed;

import basic.*;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.junit.jupiter.api.Test;
import tensorproduct.*;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.TPCell;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class MixedCellIntegralTest
{
	@Test
	public void testPressureIntegral()
	{
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(3,0),
			CoordinateVector.fromValues(5,4),new IntCoordinates(1,1));
		TPCell c = g.cells.get(0);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					ContinuousTPShapeFunction f1 = new ContinuousTPShapeFunction(c, k, i);
					ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c, k, j);
					QkQkFunction mf1 = new QkQkFunction(f1);
					QkQkFunction mf2 = new QkQkFunction(f2);
					TPCellIntegral<ContinuousTPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
						TPCellIntegral.GRAD_GRAD);
					MixedCellIntegral<TPCell,ContinuousTPShapeFunction,
						ContinuousTPVectorFunction,QkQkFunction> gg2 =
						MixedCellIntegral.fromPressureIntegral(gg);
					assertTrue(Math.abs(gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						mf1, mf2)) < 1e-12,
						"Difference of " + Math.abs(gg.evaluateCellIntegral(c, f1,
						f2) - gg2.evaluateCellIntegral(c,
						mf1, mf2)) + " for polynomial degree " + k + " and functions " + i +
							" and " + j);
				}
	}
	@Test
	public void testVelocityIntegral()
	{
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(3,0),
			CoordinateVector.fromValues(5,4), new IntCoordinates(1,1));
		TPCell c = g.cells.get(0);
		for (int k = 1; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					ContinuousTPVectorFunction f1 = new ContinuousTPVectorFunction(c, k, i);
					ContinuousTPVectorFunction f2 = new ContinuousTPVectorFunction(c, k, j);
					QkQkFunction mf1 = new QkQkFunction(f1);
					QkQkFunction mf2 = new QkQkFunction(f2);
					TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
						new TPVectorCellIntegral<>(ScalarFunction.constantFunction(1),
						TPCellIntegral.GRAD_GRAD);
					MixedCellIntegral<TPCell,ContinuousTPShapeFunction,
						ContinuousTPVectorFunction,QkQkFunction> gg2 =
						MixedCellIntegral.fromVelocityIntegral(gg);
					assertTrue(Math.abs(gg.evaluateCellIntegral(c, f1, f2) - gg2.evaluateCellIntegral(c,
						mf1, mf2)) < 1e-12,
						"Difference of " + Math.abs(gg.evaluateCellIntegral(c, f1,
							f2) - gg2.evaluateCellIntegral(c,
							mf1, mf2)) + " for polynomial degree " + k + " and functions " + i +
							" and " + j);
				}
	}
	@Test
	public void testMixedIntegral()
	{
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(3,0),
			CoordinateVector.fromValues(5,4), new IntCoordinates(1,1));
		TPCell c = g.cells.get(0);
		
		ContinuousTPVectorFunction f1 = null;
		MixedTPCellIntegral<ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vg =
			new MixedTPCellIntegral<>(MixedTPCellIntegral.VALUE_GRAD);
		MixedTPCellIntegral<ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> dv =
			new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		double [] valuesVG = new double[]{-0.666, -0.666, -0.333, -0.333, -0.333, -0.1666, -0.333, -0.1666};
		double [] valuesDV = new double[]{-0.666, 0.666, -0.333, 0.333, -0.333, -0.1666, 0.333, 0.1666};
		for(int i = 0; i < 8; i++)
		{
			System.out.println(i);
			f1 = new ContinuousTPVectorFunction(c, 1, i);
			System.out.println(f1.getNodeFunctionalPoint());
			ContinuousTPShapeFunction f2 = new ContinuousTPShapeFunction(c,1,0);
			QkQkFunction mf1 = new QkQkFunction(f1);
			QkQkFunction mf2 = new QkQkFunction(f2);
			assertTrue(Math.abs(vg.evaluateCellIntegral(c,mf1,mf2) - valuesVG[i]) < 1e-2);
			assertTrue(Math.abs(dv.evaluateCellIntegral(c,mf1,mf2) - valuesDV[i]) < 1e-2);
		}
		PlotWindow p = new PlotWindow();
		p.addPlot(new ScalarPlot2D(f1.getComponentFunction(0), g.generatePlotPoints(20),20));
		p.addPlot(new ScalarPlot2D(f1.getComponentFunction(1), g.generatePlotPoints(20),20));
		try
		{
			Thread.sleep(100000);
		} catch (InterruptedException e)
		{
			e.printStackTrace();
		}
	}
}
