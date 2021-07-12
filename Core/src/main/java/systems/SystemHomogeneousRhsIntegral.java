package systems;

import basic.Cell;
import basic.CellIntegral;
import basic.RightHandSideIntegral;
import basic.ShapeFunction;

public class SystemHomogeneousRhsIntegral<CT extends Cell<CT,?,?>, ST extends ShapeFunction<CT,?,?,?,?,?>>
	extends RightHandSideIntegral<CT, SystemShapeFunction<CT,?,?,ST>>
{
	private final RightHandSideIntegral<CT, ST> cellIntegral;
	private final int component;
	
	public SystemHomogeneousRhsIntegral(RightHandSideIntegral<CT, ST> cellIntegral, int component)
	{
		super();
		this.cellIntegral = cellIntegral;
		this.component = component;
	}
	
	
	@Override
	public double evaluateRightHandSideIntegral(CT cell, SystemShapeFunction<CT, ?, ?, ST> shapeFunction)
	{
		if(shapeFunction.mainComponent == component)
		{
			return this.cellIntegral.evaluateRightHandSideIntegral(cell,
				shapeFunction.getComponentFunction(component));
		}
		return 0;
	}
}
