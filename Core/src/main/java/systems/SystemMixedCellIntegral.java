package systems;

import basic.*;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.QuadratureRule1D;
import tensorproduct.TPCell;
import tensorproduct.TPCellIntegral;

public class SystemMixedCellIntegral<CT extends Cell<CT,?,?>>
	extends CellIntegral<CT, SystemShapeFunction<CT,?,?,?>>
{
	private final int component1;
	private final int component2;
	public static final String VALUE_VALUE = "ValueValue";
	public static final String GRAD_GRAD = "GradGrad";
	public final QuadratureRule1D quadratureRule1D;
	
	public SystemMixedCellIntegral(String name, int component1, int component2)
	{
		this(name, component1, component2, QuadratureRule1D.Gauss5);
	}
	
	public SystemMixedCellIntegral(String name, int component1, int component2, QuadratureRule1D quadratureRule1D)
	{
		super(name, quadratureRule1D);
		this.component1 = component1;
		this.component2 = component2;
		this.quadratureRule1D = quadratureRule1D;
	}
	
	public SystemMixedCellIntegral(String name, Function<?, ?, ?> weight, int component1, int component2)
	{
		this(name, weight, component1, component2, QuadratureRule1D.Gauss5);
	}
	
	public SystemMixedCellIntegral(String name, Function<?, ?, ?> weight, int component1, int component2, QuadratureRule1D quadratureRule1D)
	{
		super(weight, name, quadratureRule1D);
		this.component1 = component1;
		this.component2 = component2;
		this.quadratureRule1D = quadratureRule1D;
	}
	
	@Override
	public double evaluateCellIntegral(CT cell, SystemShapeFunction<CT, ?, ?, ?> shapeFunction1, SystemShapeFunction<CT, ?, ?, ?> shapeFunction2)
	{
		if (shapeFunction1.mainComponent == component1 && shapeFunction2.mainComponent == component2)
		{
			if (name.equals(VALUE_VALUE))
			{
				if (Double.class.isAssignableFrom(SystemParameters.getInstance().signatures[component1].getValueT()))
					return TPCellIntegral.integrateNonTensorProduct(x ->
							(Double) shapeFunction1.getComponentFunction(component1).value(x)
								* (Double) shapeFunction2.getComponentFunction(component2).value(x) * (Double) weight.value(x),
						((TPCell) cell).getCell1Ds(), quadratureRule1D);
				if (CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[component1].getValueT()))
					return TPCellIntegral.integrateNonTensorProduct(x ->
							((CoordinateVector) shapeFunction1.getComponentFunction(component1).value(x))
								.inner((CoordinateVector) shapeFunction2.getComponentFunction(component2).value(x)) * (Double) weight.value(x),
						((TPCell) cell).getCell1Ds(), quadratureRule1D);
			}
			if (name.equals(GRAD_GRAD))
			{
				if (CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[component1].getGradientT()))
					return TPCellIntegral.integrateNonTensorProduct(x ->
							((CoordinateVector) shapeFunction1.getComponentFunction(component1).gradient(x))
								.inner((CoordinateVector) shapeFunction2.getComponentFunction(component2).gradient(x)) * (Double) weight.value(x),
						((TPCell) cell).getCell1Ds(), quadratureRule1D);
				if (CoordinateMatrix.class.isAssignableFrom(SystemParameters.getInstance().signatures[component1].getGradientT()))
					return TPCellIntegral.integrateNonTensorProduct(x ->
							((CoordinateMatrix) shapeFunction1.getComponentFunction(component1).gradient(x))
								.frobeniusInner((CoordinateMatrix) shapeFunction2.getComponentFunction(component2).gradient(x)) * (Double) weight.value(x),
						((TPCell) cell).getCell1Ds(), quadratureRule1D);
			}
		}
		return 0;
	}
	
}
