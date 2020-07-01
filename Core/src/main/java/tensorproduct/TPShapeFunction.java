package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.ScalarShapeFunction;
import basic.Cell;
import linalg.DoubleTensor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class TPShapeFunction extends ScalarShapeFunction
{
	LagrangeBasisFunction1D xFunction;
	LagrangeBasisFunction1D yFunction;
	private DoubleTensor functional_point;
	public TPShapeFunction(LagrangeBasisFunction1D xFunction, LagrangeBasisFunction1D yFunction,TPCell cell)
	{
		this.xFunction = xFunction;
		this.yFunction = yFunction;
		functional_point = DoubleTensor.vectorFromValues(xFunction.degreeOfFreedom, yFunction.degreeOfFreedom);
		this.cells = new ArrayList<>();
		cells.add(cell);
		nodeFunctional = new LagrangeNodeFunctional(functional_point);
	}

	@Override
	public double value(DoubleTensor pos)
	{
		return xFunction.value(pos.x())*yFunction.value(pos.y());
	}

	@Override
	public DoubleTensor derivative(DoubleTensor pos)
	{
		return DoubleTensor.vectorFromValues(xFunction.derivative(pos.x())*yFunction.value(pos.y()),
			xFunction.value(pos.x())*yFunction.derivative(pos.y()));
	}

	@Override
	public double valueInCell(DoubleTensor pos, Cell cell)
	{
		if(cells.get(0) == cell)
			return value(pos);
		else
			return 0;
	}
	@Override
	public DoubleTensor derivativeInCell(DoubleTensor pos, Cell cell)
	{
		if(cells.get(0) == cell)
			return derivative(pos);
		else
			return new DoubleTensor(pos.size());
	}

	@Override
	public Map<Integer, Double> prolongate(ArrayList<ScalarShapeFunction> refinedFunctions)
	{
		Map<Integer, Double> ret = new HashMap<>();
		for(ScalarShapeFunction shapeFunction:refinedFunctions)
		{
			ret.put(shapeFunction.getGlobalIndex(), shapeFunction.getNodeFunctional().evaluate(this));
		}
		return ret;
	}
}
