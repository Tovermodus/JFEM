package tensorproduct;

import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.*;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

public class IntegralTest
{
	
	@Test
	public void testTPCellIntegralDeg1G3() {
		int polynomialDegree = 1;
		double[] testValuesgg = new double[]{0.683,-0.283,-0.058,-0.342,-0.283,0.683,-0.342,-0.058,-0.058,-0.342,0.683,-0.283,-0.342,-0.058,-0.283,0.683};
		double[] testValuesggw = new double[]{1.564, -0.251, -0.168, -1.144, -0.251, 1.617, -1.144, -0.222, -0.168, -1.144, 2.96, -1.647, -1.144, -0.222, -1.647, 3.013};
		double[] testValuesvv = new double[]{0.088, 0.044, 0.044, 0.022, 0.044, 0.088, 0.022, 0.044, 0.044, 0.022, 0.088, 0.044, 0.022, 0.044, 0.044, 0.088};
		double[] testValuesvvw = new double[]{0.144, 0.074, 0.144, 0.074, 0.074, 0.153, 0.074, 0.153, 0.144, 0.074, 0.433, 0.223, 0.074, 0.153, 0.223, 0.46};
		double[] testValuesvg = new double[]{-0.02, -0.15, 0.154, 0.016, 0.074, -0.263, 0.125, 0.063, -0.179, -0.163, 0.304, 0.038, -0.054, -0.323, 0.246, 0.131};
		double[] testValuesgv = new double[]{-0.02, 0.074, -0.179, -0.054, -0.15, -0.263, -0.163, -0.323, 0.154, 0.125, 0.304, 0.246, 0.016, 0.063, 0.038, 0.131};
		
		Cell1D c2 = new Cell1D(0,1,QuadratureRule1D.Gauss3);
		Cell1D c1 = new Cell1D(2.3,3.1,QuadratureRule1D.Gauss3);
		TPCell cell = new TPCell(List.of(c1,c2));
		List<TPShapeFunction> shapeFunctions = new ArrayList<>();
		for (int i = 0; i < Math.pow(polynomialDegree+1,2); i++)
		{
			shapeFunctions.add(new TPShapeFunction(cell, polynomialDegree, i));
		}
		ScalarFunction weight = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return pos.at(1)*(pos.at(0)+4);
			}
		};
		VectorFunction vweight = new VectorFunction()
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
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(2-pos.x()+pos.y()*0.1,pos.x()*pos.y());
			}
		};
		TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		TPCellIntegral<TPShapeFunction> ggw = new TPCellIntegral<>(weight,TPCellIntegral.GRAD_GRAD);
		TPCellIntegral<TPShapeFunction> vg = new TPCellIntegral<>(vweight,TPCellIntegral.VALUE_GRAD);
		TPCellIntegral<TPShapeFunction> gv = new TPCellIntegral<>(vweight,TPCellIntegral.GRAD_VALUE);
		TPCellIntegral<TPShapeFunction> vv = new TPCellIntegral<>(TPCellIntegral.VALUE_VALUE);
		TPCellIntegral<TPShapeFunction> vvw = new TPCellIntegral<>(weight,TPCellIntegral.VALUE_VALUE);
		int k = 0;
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				assertTrue(Math.abs(gg.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesgg[k])<=1e-2);
				assertTrue(Math.abs(ggw.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesggw[k])<=1e-2);
				assertTrue(Math.abs(vg.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesvg[k])<=1e-2);
				assertTrue(Math.abs(gv.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesgv[k])<=1e-2);
				assertTrue(Math.abs(vv.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesvv[k])<=1e-2);
				assertTrue(Math.abs(vvw.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesvvw[k++])<=1e-2);
			}
		}
		
	}
	
	@Test
	public void testTPCellIntegralDeg1G5() {
		int polynomialDegree = 1;
		double[] testValuesgg = new double[]{0.683,-0.283,-0.058,-0.342,-0.283,0.683,-0.342,-0.058,-0.058,-0.342,0.683,-0.283,-0.342,-0.058,-0.283,0.683};
		double[] testValuesggw = new double[]{1.564, -0.251, -0.168, -1.144, -0.251, 1.617, -1.144, -0.222, -0.168, -1.144, 2.96, -1.647, -1.144, -0.222, -1.647, 3.013};
		double[] testValuesvv = new double[]{0.088, 0.044, 0.044, 0.022, 0.044, 0.088, 0.022, 0.044, 0.044, 0.022, 0.088, 0.044, 0.022, 0.044, 0.044, 0.088};
		double[] testValuesvvw = new double[]{0.144, 0.074, 0.144, 0.074, 0.074, 0.153, 0.074, 0.153, 0.144, 0.074, 0.433, 0.223, 0.074, 0.153, 0.223, 0.46};
		double[] testValuesvg = new double[]{-0.02, -0.15, 0.154, 0.016, 0.074, -0.263, 0.125, 0.063, -0.179, -0.163, 0.304, 0.038, -0.054, -0.323, 0.246, 0.131};
		double[] testValuesgv = new double[]{-0.02, 0.074, -0.179, -0.054, -0.15, -0.263, -0.163, -0.323, 0.154, 0.125, 0.304, 0.246, 0.016, 0.063, 0.038, 0.131};
		
		Cell1D c2 = new Cell1D(0,1,QuadratureRule1D.Gauss5);
		Cell1D c1 = new Cell1D(2.3,3.1,QuadratureRule1D.Gauss5);
		TPCell cell = new TPCell(List.of(c1,c2));
		List<TPShapeFunction> shapeFunctions = new ArrayList<>();
		for (int i = 0; i < Math.pow(polynomialDegree+1,2); i++)
		{
			shapeFunctions.add(new TPShapeFunction(cell, polynomialDegree, i));
		}
		ScalarFunction weight = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return pos.at(1)*(pos.at(0)+4);
			}
		};
		VectorFunction vweight = new VectorFunction()
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
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(2-pos.x()+pos.y()*0.1,pos.x()*pos.y());
			}
		};
		TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		TPCellIntegral<TPShapeFunction> ggw = new TPCellIntegral<>(weight,TPCellIntegral.GRAD_GRAD);
		TPCellIntegral<TPShapeFunction> vg = new TPCellIntegral<>(vweight,TPCellIntegral.VALUE_GRAD);
		TPCellIntegral<TPShapeFunction> gv = new TPCellIntegral<>(vweight,TPCellIntegral.GRAD_VALUE);
		TPCellIntegral<TPShapeFunction> vv = new TPCellIntegral<>(TPCellIntegral.VALUE_VALUE);
		TPCellIntegral<TPShapeFunction> vvw = new TPCellIntegral<>(weight,TPCellIntegral.VALUE_VALUE);
		int k = 0;
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				assertTrue(Math.abs(gg.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesgg[k])<=1e-2);
				assertTrue(Math.abs(ggw.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesggw[k])<=1e-2);
				assertTrue(Math.abs(vg.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesvg[k])<=1e-2);
				assertTrue(Math.abs(gv.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesgv[k])<=1e-2);
				assertTrue(Math.abs(vv.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesvv[k])<=1e-2);
				assertTrue(Math.abs(vvw.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesvvw[k++])<=1e-2);
			}
		}
		
	}
	
	@Test
	public void testTPCellIntegralDeg2G5() {
		int polynomialDegree = 2;
		double[] testValuesgg = new double[]{0.637, -0.319, -0.006, -0.089, -0.364, 0.098, -0.061, 0.128, -0.022, -0.319, 1.884, -0.32, -0.364, -0.693, -0.364, 0.128, -0.08, 0.128, -0.006, -0.32, 0.637, 0.098, -0.364, -0.09, -0.022, 0.128, -0.061, -0.089, -0.364, 0.098, 2.124, -1.493, 0.079, -0.089, -0.364, 0.098, -0.364, -0.693, -0.364, -1.493, 5.831, -1.493, -0.364, -0.693, -0.364, 0.098, -0.364, -0.09, 0.079, -1.493, 2.124, 0.098, -0.364, -0.09, -0.061, 0.128, -0.022, -0.089, -0.364, 0.098, 0.637, -0.319, -0.006, 0.128, -0.08, 0.128, -0.364, -0.693, -0.364, -0.319, 1.884, -0.319, -0.022, 0.128, -0.061, 0.098, -0.364, -0.09, -0.006, -0.319, 0.637};
		double[] testValuesggw = new double[]{0.655, -0.193, -0.042, -0.455, -0.223, 0.119, -0.2, 0.417, -0.076, -0.193, 2.173, -0.194, -0.223, -1.905, -0.252, 0.417, -0.268, 0.446, -0.042, -0.194, 0.71, 0.119, -0.252, -0.497, -0.076, 0.446, -0.212, -0.455, -0.223, 0.119, 6.853, -4.881, 0.267, -0.106, -2.116, 0.543, -0.223, -1.905, -0.252, -4.881, 19.534, -5.123, -2.116, -2.739, -2.29, 0.119, -0.252, -0.497, 0.267, -5.123, 7.38, 0.543, -2.29, -0.146, -0.2, 0.417, -0.076, -0.106, -2.116, 0.543, 3.453, -1.911, -0.001, 0.417, -0.268, 0.446, -2.116, -2.739, -2.29, -1.911, 10.451, -1.989, -0.076, 0.446, -0.212, 0.543, -2.29, -0.146, -0.001, -1.989, 3.726};
		double[] testValuesvv = new double[]{0.014, 0.007, -0.003, 0.007, 0.003, -0.001, -0.003, -0.001, 0.0, 0.007, 0.056, 0.007, 0.003, 0.028, 0.003, -0.001, -0.014, -0.001, -0.003, 0.007, 0.014, -0.001, 0.003, 0.007, 0.0, -0.001, -0.003, 0.007, 0.003, -0.001, 0.056, 0.028, -0.014, 0.007, 0.003, -0.001, 0.003, 0.028, 0.003, 0.028, 0.227, 0.028, 0.003, 0.028, 0.003, -0.001, 0.003, 0.007, -0.014, 0.028, 0.056, -0.001, 0.003, 0.007, -0.003, -0.001, 0.0, 0.007, 0.003, -0.001, 0.014, 0.007, -0.003, -0.001, -0.014, -0.001, 0.003, 0.028, 0.003, 0.007, 0.056, 0.007, 0.0, -0.001, -0.003, -0.001, 0.003, 0.007, -0.003, 0.007, 0.014};
		double[] testValuesvvw = new double[]{0.011, 0.005, -0.002, 0.0, 0.0, 0.0, -0.011, -0.005, 0.002, 0.005, 0.047, 0.006, 0.0, 0.0, 0.0, -0.005, -0.047, -0.006, -0.002, 0.006, 0.012, 0.0, 0.0, 0.0, 0.002, -0.006, -0.012, 0.0, 0.0, 0.0, 0.182, 0.089, -0.047, 0.045, 0.022, -0.011, 0.0, 0.0, 0.0, 0.089, 0.762, 0.1, 0.022, 0.19, 0.025, 0.0, 0.0, 0.0, -0.047, 0.1, 0.199, -0.011, 0.025, 0.049, -0.011, -0.005, 0.002, 0.045, 0.022, -0.011, 0.079, 0.039, -0.02, -0.005, -0.047, -0.006, 0.022, 0.19, 0.025, 0.039, 0.333, 0.044, 0.002, -0.006, -0.012, -0.011, 0.025, 0.049, -0.02, 0.044, 0.087};
		double[] testValuesvg = new double[]{0.009, -0.047, 0.018, 0.047, -0.004, -0.002, -0.023, 0.0, 0.001, 0.038, -0.048, -0.086, 0.04, 0.167, -0.016, -0.019, -0.083, 0.006, -0.012, 0.071, -0.086, -0.018, 0.063, 0.009, 0.008, -0.03, -0.005, -0.037, -0.044, 0.021, 0.026, -0.178, 0.069, 0.129, 0.041, -0.027, 0.0, -0.216, -0.071, 0.141, -0.193, -0.332, 0.076, 0.551, 0.043, 0.005, 0.008, -0.097, -0.045, 0.272, -0.336, -0.041, 0.114, 0.119, 0.019, 0.021, -0.01, -0.126, -0.081, 0.044, 0.132, 0.02, -0.019, 0.001, 0.108, 0.034, -0.045, -0.6, -0.121, 0.093, 0.527, 0.002, -0.003, -0.003, 0.047, 0.03, -0.05, -0.2, -0.046, 0.147, 0.078};
		double[] testValuesgv = new double[]{0.009, 0.038, -0.012, -0.037, 0.0, 0.005, 0.019, 0.001, -0.003, -0.047, -0.048, 0.071, -0.044, -0.216, 0.008, 0.021, 0.108, -0.003, 0.018, -0.086, -0.086, 0.021, -0.071, -0.097, -0.01, 0.034, 0.047, 0.047, 0.04, -0.018, 0.026, 0.141, -0.045, -0.126, -0.045, 0.03, -0.004, 0.167, 0.063, -0.178, -0.193, 0.272, -0.081, -0.6, -0.05, -0.002, -0.016, 0.009, 0.069, -0.332, -0.336, 0.044, -0.121, -0.2, -0.023, -0.019, 0.008, 0.129, 0.076, -0.041, 0.132, 0.093, -0.046, 0.0, -0.083, -0.03, 0.041, 0.551, 0.114, 0.02, 0.527, 0.147, 0.001, 0.006, -0.005, -0.027, 0.043, 0.119, -0.019, 0.002, 0.078};
		
		Cell1D c2 = new Cell1D(0,1,QuadratureRule1D.Gauss5);
		Cell1D c1 = new Cell1D(2.3,3.1,QuadratureRule1D.Gauss5);
		TPCell cell = new TPCell(List.of(c1,c2));
		List<TPShapeFunction> shapeFunctions = new ArrayList<>();
		for (int i = 0; i < Math.pow(polynomialDegree+1,2); i++)
		{
			shapeFunctions.add(new TPShapeFunction(cell, polynomialDegree, i));
		}
		ScalarFunction weight = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return pos.at(1)*(pos.at(0)+4);
			}
		};
		VectorFunction vweight = new VectorFunction()
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
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(2-pos.x()+pos.y()*0.1,pos.x()*pos.y());
			}
		};
		TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		TPCellIntegral<TPShapeFunction> ggw = new TPCellIntegral<>(weight,TPCellIntegral.GRAD_GRAD);
		TPCellIntegral<TPShapeFunction> vg = new TPCellIntegral<>(vweight,TPCellIntegral.VALUE_GRAD);
		TPCellIntegral<TPShapeFunction> gv = new TPCellIntegral<>(vweight,TPCellIntegral.GRAD_VALUE);
		TPCellIntegral<TPShapeFunction> vv = new TPCellIntegral<>(TPCellIntegral.VALUE_VALUE);
		TPCellIntegral<TPShapeFunction> vvw = new TPCellIntegral<>(weight,TPCellIntegral.VALUE_VALUE);
		int k = 0;
		for (int i = 0; i < shapeFunctions.size(); i++)
		{
			for (int j = 0; j < shapeFunctions.size(); j++)
			{
				assertTrue(Math.abs(gg.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesgg[k])<=1e-2);
				assertTrue(Math.abs(ggw.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesggw[k])<=1e-2);
				assertTrue(Math.abs(vg.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesvg[k])<=1e-2);
				assertTrue(Math.abs(gv.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesgv[k])<=1e-2);
				assertTrue(Math.abs(vv.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesvv[k])<=1e-2);
				assertTrue(Math.abs(vvw.evaluateCellIntegral(cell, shapeFunctions.get(i),
					shapeFunctions.get(j))-testValuesvvw[k++])<=1e-2);
			}
		}
		
	}
	
	
}