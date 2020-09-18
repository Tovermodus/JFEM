package tensorproduct;

import basic.CellIntegral;
import basic.Function;
import basic.ScalarShapeFunction;
import com.google.common.primitives.Doubles;
import linalg.CoordinateVector;
import linalg.Vector;

import java.util.ArrayList;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class TPCellIntegral<ST extends ScalarShapeFunction<TPCell,TPFace,ST>> extends CellIntegral<TPCell,
	TPFace
	,ST>
{
	public static final String GRAD_GRAD = "GradGrad";
	public static final String VALUE_VALUE = "ValueValue";
	public static final String GRAD_VALUE = "GradValue";
	public static final String VALUE_GRAD = "ValueGrad";
	private final boolean weightIsTensorProduct;
	public TPCellIntegral(Function<?,?,?> weight, String name, boolean weightIsTensorProduct)
	{
		super(weight,name);
		this.weightIsTensorProduct = weightIsTensorProduct;
		if(name.equals(GRAD_GRAD) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(GRAD_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Vector))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_GRAD) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Vector))
			throw new IllegalArgumentException();
			
	}
	public TPCellIntegral(String name)
	{
		super(name);
		this.weightIsTensorProduct = true;
		if(name.equals(GRAD_GRAD) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(GRAD_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Vector))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_GRAD) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Vector))
			throw new IllegalArgumentException();
	}
	public static double integrateNonTensorProduct(ToDoubleFunction<CoordinateVector> eval, List<Cell1D> cells)
	{
		int qsize = cells.get(0).points.length;
		double ret = 0;
		double val = 0;
		CoordinateVector quadraturePoint = new CoordinateVector(cells.size());
		int [] decomposedPointIndex = new int[cells.size()];
		for (int i = 0; i < Math.pow(qsize,cells.size()); i++)
		{
			int icopy = i;
			for (int j = 0; j < cells.size(); j++)
			{
				decomposedPointIndex[j] = icopy % qsize;
				icopy = icopy/qsize;
				quadraturePoint.set(cells.get(j).points[decomposedPointIndex[j]],j);
			}
			//System.out.println(quadraturePoint);
			val = eval.applyAsDouble(quadraturePoint);
			//System.out.println(ret + " ret-: " + val+ " "+ quadraturePoint);
			for(int j = 0; j < cells.size(); j++)
				val *= cells.get(j).weights[decomposedPointIndex[j]];
			ret += val;
		}
		return ret;
	}
	@Override
	public double evaluateCellIntegral(TPCell cell, ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		if(name.equals(GRAD_GRAD))
		{
				return integrateNonTensorProduct(x->shapeFunction1.gradient(x).inner(shapeFunction2.gradient(x))*(Double)weight.value(x),cell.cell1Ds);
		}
		if(name.equals(VALUE_VALUE))
		{
				return integrateNonTensorProduct(x->shapeFunction1.value(x)*shapeFunction2.value(x)*(Double)weight.value(x),cell.cell1Ds);
		}
		if(name.equals(GRAD_VALUE))
		{
				return integrateNonTensorProduct(x->shapeFunction1.gradient(x).inner((Vector)weight.value(x))*shapeFunction2.value(x),cell.cell1Ds);
		}
		if(name.equals(VALUE_GRAD))
		{
				return integrateNonTensorProduct(x->shapeFunction2.gradient(x).inner((Vector)weight.value(x))*shapeFunction1.value(x),cell.cell1Ds);
		}
		throw new UnsupportedOperationException("unknown integral name");
	}
}