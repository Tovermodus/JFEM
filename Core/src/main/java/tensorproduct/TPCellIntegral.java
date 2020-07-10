package tensorproduct;

import basic.CellIntegral;
import basic.Function;
import linalg.CoordinateVector;
import linalg.Vector;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class TPCellIntegral extends CellIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>
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
	static double integrateTensorProduct(ToDoubleFunction<CoordinateVector> eval, List<Cell1D> cells)
	{
		double ret = 1;
		for(int j = 0; j < cells.size(); j++)
		{
			double val = 0;
			Cell1D cell = cells.get(j);
			CoordinateVector quadPoint =
				CoordinateVector.fromValues(cells.stream().mapToDouble(Cell1D::center).toArray());
			for(int i = 0; i < cell.points.length; i++)
			{
				quadPoint.set(cell.points[i],j);
				val+= eval.applyAsDouble(quadPoint)*cell.weights[i];
			}
			ret *= val;
		}
		return ret;
	}
	static double integrateNonTensorProduct(ToDoubleFunction<CoordinateVector> eval, List<Cell1D> cells)
	{
		if(cells.size() == 1)
		{
			double ret = 0;
			Cell1D cell = cells.get(0);
			for(int i = 0; i < cell.points.length; i++)
			{
				ret += eval.applyAsDouble(CoordinateVector.fromValues(cell.points[i]))*cell.weights[i];
			}
			return ret;
		}
		
		double ret = 0;
		Cell1D cell = cells.get(0);
		for(int i = 0; i < cell.points.length; i++)
		{
			int finalI = i;
			ret += integrateNonTensorProduct(x->{
				double [] point = new double [cells.size()];
				point[0] = cell.points[finalI];
				for (int j = 1; j < cells.size(); j++)
				{
					point[j] = x.at(j-1);
				}
				return eval.applyAsDouble(CoordinateVector.fromValues(point));
			},cells.subList(1,cells.size()-1))*cell.weights[i];
		}
		return ret;
	}
	@Override
	public double evaluateCellIntegral(TPCell<TPShapeFunction> cell, TPShapeFunction shapeFunction1, TPShapeFunction shapeFunction2)
	{
		if(name.equals(GRAD_GRAD))
		{
			if(weightIsTensorProduct)
				return integrateTensorProduct(x->shapeFunction1.gradient(x).inner(shapeFunction2.gradient(x))*(Double)weight.value(x),cell.cell1Ds);
			else
				return integrateNonTensorProduct(x->shapeFunction1.gradient(x).inner(shapeFunction2.gradient(x))*(Double)weight.value(x),cell.cell1Ds);
		}
		if(name.equals(VALUE_VALUE))
		{
			if(weightIsTensorProduct)
				return integrateTensorProduct(x->shapeFunction1.value(x)*shapeFunction2.value(x)*(Double)weight.value(x),cell.cell1Ds);
			else
				return integrateNonTensorProduct(x->shapeFunction1.value(x)*shapeFunction2.value(x)*(Double)weight.value(x),cell.cell1Ds);
		}
		if(name.equals(GRAD_VALUE))
		{
			if(weightIsTensorProduct)
				return integrateTensorProduct(x->shapeFunction1.gradient(x).inner((Vector)weight.value(x))*shapeFunction2.value(x),cell.cell1Ds);
			else
				return integrateNonTensorProduct(x->shapeFunction1.gradient(x).inner((Vector)weight.value(x))*shapeFunction2.value(x),cell.cell1Ds);
		}
		if(name.equals(VALUE_GRAD))
		{
			if(weightIsTensorProduct)
				return integrateTensorProduct(x->shapeFunction2.gradient(x).inner((Vector)weight.value(x))*shapeFunction1.value(x),cell.cell1Ds);
			else
				return integrateNonTensorProduct(x->shapeFunction2.gradient(x).inner((Vector)weight.value(x))*shapeFunction1.value(x),cell.cell1Ds);
		}
		throw new UnsupportedOperationException("unknown integral name");
	}
}