package systems;

import basic.*;

public class SystemRightHandSideIntegral<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>>
	extends RightHandSideIntegral<CT, SystemShapeFunction<CT, FT, ET>>
{
	@Override
	public double evaluateRightHandSideIntegral(CT cell, SystemShapeFunction<CT, FT, ET> shapeFunction1)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
