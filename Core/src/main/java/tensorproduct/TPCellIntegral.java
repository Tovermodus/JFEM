package tensorproduct;

import basic.CellIntegral;
import basic.ScalarFunction;
import basic.ScalarShapeFunction;
import basic.TensorFunction;
import basic.Cell;
import linalg.DoubleTensor;
import tensorproduct.TPFunction;

public class TPCellIntegral extends CellIntegral
{
	TPFunction tpWeight;
	TPVectorFunction tpVectorWeight;
	public static String GRAD_GRAD = "GradGrad";
	public static String VALUE_VALUE = "ValueValue";
	public static String GRAD_VALUE = "GradValue";
	public static String VALUE_GRAD = "ValueGrad";

	public TPCellIntegral(ScalarFunction vectorWeight, String name)
	{
		super(vectorWeight, name);
	}
	public TPCellIntegral(TensorFunction weight, String name)
	{
		super(weight, name);
	}

	public TPCellIntegral(TPFunction weight, String name)
	{
		super(weight, name);
		this.tpWeight = weight;
	}


	public TPCellIntegral(TPVectorFunction vectorWeight, String name)
	{
		super(vectorWeight, name);
		this.tpVectorWeight = vectorWeight;
	}

	public TPCellIntegral(String name)
	{
		this(new TPFunction(Function1D.oneFunction(),Function1D.oneFunction()),name);
	}

	@Override
	public double evaluateCellIntegral(Cell cell, ScalarShapeFunction shapeFunction1, ScalarShapeFunction shapeFunction2)
	{
		if(cell instanceof TPCell)
		{
			if(shapeFunction1 instanceof TPShapeFunction && shapeFunction2 instanceof TPShapeFunction && (tpWeight != null || tpVectorWeight != null))
				return this.evaluateCellIntegral((TPCell) cell,(TPShapeFunction)shapeFunction1,
					(TPShapeFunction)shapeFunction2);
			else
				return this.evaluateCellIntegral((TPCell) cell,shapeFunction1,
					shapeFunction2);

		}
		throw new UnsupportedOperationException();
	}

	public double evaluateCellIntegral(TPCell cell, ScalarShapeFunction function1,
	                                   ScalarShapeFunction function2)
	{
		double ret = 0;
		for(int i = 0; i < cell.cellx.weights.length; i++)
			for(int j = 0; j < cell.celly.weights.length; j++)
			{
				DoubleTensor point = DoubleTensor.vectorFromValues(cell.cellx.points[i],
					cell.celly.points[j]);
				double integralWeight = cell.cellx.weights[i]*cell.celly.weights[j];
				if(name.equals(GRAD_GRAD))
					ret += function1.derivative(point).inner(function2.derivative(point))*weight.value(point)*integralWeight;
				if(name.equals(VALUE_VALUE))
					ret += function1.value(point)*function2.value(point)*weight.value(point)*integralWeight;
				if(name.equals(GRAD_VALUE))
					ret += function1.derivative(point).inner(vectorWeight.value(point))*function2.value(point)*integralWeight;
				if(name.equals(VALUE_GRAD))
					ret += function1.value(point)*function2.derivative(point).inner(vectorWeight.value(point))*integralWeight;
			}
		return ret;
	}
	double evaluateCellIntegral(TPCell cell, TPShapeFunction function1,
	                             TPShapeFunction function2)
	{
		Cell1D cellx = cell.cellx;
		Cell1D celly = cell.celly;
		if(name.equals(GRAD_GRAD))
		{
			double derX = 0;
			double derY = 0;
			double valX = 0;
			double valY = 0;
			for(int i = 0; i < cellx.weights.length; i++)
			{
				double αi = cellx.weights[i];
				double xi = cellx.points[i];
				derX += function1.xFunction.derivative(xi)*function2.xFunction.derivative(xi)*this.tpWeight.xFunction.value(xi)*αi;
				valX += function1.xFunction.value(xi)*function2.xFunction.value(xi)*this.tpWeight.xFunction.value(xi)*αi;
			}
			for(int i = 0; i < celly.weights.length; i++)
			{
				double αi = celly.weights[i];
				double yi = celly.points[i];
				derY += function1.yFunction.derivative(yi)*function2.yFunction.derivative(yi)*this.tpWeight.yFunction.value(yi)*αi;
				valY += function1.yFunction.value(yi)*function2.yFunction.value(yi)*this.tpWeight.yFunction.value(yi)*αi;
			}
			return derX*valY+valX*derY;
		}
		if(name.equals(VALUE_VALUE))
		{
			double valX = 0;
			double valY = 0;
			for(int i = 0; i < cellx.weights.length; i++)
			{
				double αi = cellx.weights[i];
				double xi = cellx.points[i];
				valX += function1.xFunction.value(xi)*function2.xFunction.value(xi)*this.tpWeight.xFunction.value(xi)*αi;
			}
			for(int i = 0; i < celly.weights.length; i++)
			{
				double αi = celly.weights[i];
				double yi = celly.points[i];
				valY += function1.yFunction.value(yi)*function2.yFunction.value(yi)*this.tpWeight.yFunction.value(yi)*αi;
			}
			return valX*valY;
		}
		if(name.equals(GRAD_VALUE))
		{
			double derX = 0;
			double derY = 0;
			double valX = 0;
			double valY = 0;
			for(int i = 0; i < cellx.weights.length; i++)
			{
				double αi = cellx.weights[i];
				double xi = cellx.points[i];
				derX += function1.xFunction.derivative(xi)*function2.xFunction.value(xi)*this.tpVectorWeight.xFunctionx.value(xi)*αi;
				valX += function1.xFunction.value(xi)*function2.xFunction.value(xi)*this.tpVectorWeight.yFunctionx.value(xi)*αi;
			}
			for(int i = 0; i < celly.weights.length; i++)
			{
				double αi = celly.weights[i];
				double yi = celly.points[i];
				derY += function1.yFunction.derivative(yi)*function2.yFunction.value(yi)*this.tpVectorWeight.yFunctiony.value(yi)*αi;
				valY += function1.yFunction.value(yi)*function2.yFunction.value(yi)*this.tpVectorWeight.xFunctiony.value(yi)*αi;
			}
			return derX*valY+valX*valY;
		}
		if(name.equals(VALUE_GRAD))
		{
			double derX = 0;
			double derY = 0;
			double valX = 0;
			double valY = 0;
			for(int i = 0; i < cellx.weights.length; i++)
			{
				double αi = cellx.weights[i];
				double xi = cellx.points[i];
				derX += function1.xFunction.value(xi)*function2.xFunction.derivative(xi)*this.tpVectorWeight.xFunctionx.value(xi)*αi;
				valX += function1.xFunction.value(xi)*function2.xFunction.value(xi)*this.tpVectorWeight.yFunctionx.value(xi)*αi;
			}
			for(int i = 0; i < celly.weights.length; i++)
			{
				double αi = celly.weights[i];
				double yi = celly.points[i];
				derY += function1.yFunction.value(yi)*function2.yFunction.derivative(yi)*this.tpVectorWeight.yFunctiony.value(yi)*αi;
				valY += function1.yFunction.value(yi)*function2.yFunction.value(yi)*this.tpVectorWeight.xFunctiony.value(yi)*αi;
			}
			return derX*valY+valX*valY;
		}
		return 0;
	}
}