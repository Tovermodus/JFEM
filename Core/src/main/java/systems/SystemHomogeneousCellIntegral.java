package systems;

import basic.Cell;
import basic.CellIntegral;
import basic.ShapeFunction;

public class SystemHomogeneousCellIntegral<CT extends Cell<CT,?,?>, ST extends ShapeFunction<CT,?,?,?,?,?>>
	extends CellIntegral<CT, SystemShapeFunction<CT,?,?,ST>>
{
	private final CellIntegral<CT, ST> cellIntegral;
	private final int component;
	
	public SystemHomogeneousCellIntegral(CellIntegral<CT, ST> cellIntegral, int component)
	{
		super("homogeneous");
		this.cellIntegral = cellIntegral;
		this.component = component;
	}
	
	
	@Override
	public double evaluateCellIntegral(CT cell, SystemShapeFunction<CT, ?, ?, ST> shapeFunction1,
	                                   SystemShapeFunction<CT, ?, ?, ST> shapeFunction2)
	{
		if(shapeFunction1.mainComponent == component && shapeFunction2.mainComponent == component)
		{
			return this.cellIntegral.evaluateCellIntegral(cell,
				shapeFunction1.getComponentFunction(component),
				shapeFunction2.getComponentFunction(component));
		}
		return 0;
	}
	
}
