package systems;

import basic.*;

public class SystemCellIntegral<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>>
	extends CellIntegral<CT, SystemShapeFunction<CT, FT, ET>>
{
	@Override
	public double evaluateCellIntegral(CT cell, SystemShapeFunction<CT, FT, ET> shapeFunction1, SystemShapeFunction<CT, FT, ET> shapeFunction2)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
