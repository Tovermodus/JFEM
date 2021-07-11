package systems;

import basic.*;
import tensorproduct.TPCell;
import tensorproduct.TPCellIntegral;

public class SystemMixedCellIntegral<CT extends Cell<CT,?,?>>
	extends CellIntegral<CT, SystemShapeFunction<CT,?,?,?>>
{
	public SystemMixedCellIntegral(String name, int component1, int component2)
	{
	
	}
	public SystemMixedCellIntegral(){}
	
	@Override
	public double evaluateCellIntegral(CT cell, SystemShapeFunction<CT, ?, ?, ?> shapeFunction1, SystemShapeFunction<CT, ?, ?, ?> shapeFunction2)
	{
		return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x).getComponent(0).at(0)*shapeFunction2.value(x).getComponent(3).euclidianNorm(),
			((TPCell) cell).getCell1Ds());
	}
	
}
