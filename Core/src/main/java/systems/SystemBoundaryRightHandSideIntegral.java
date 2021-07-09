package systems;

import basic.*;

public class SystemBoundaryRightHandSideIntegral<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>>
	extends BoundaryRightHandSideIntegral<FT, SystemShapeFunction<CT, FT, ET>>
{
	@Override
	public double evaluateBoundaryRightHandSideIntegral(FT face, SystemShapeFunction<CT, FT, ET> shapeFunction1)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
