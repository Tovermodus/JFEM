package systems;

import basic.*;

public class SystemFaceIntegral <CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>>
	extends FaceIntegral<FT, SystemShapeFunction<CT, FT, ET>>
{
	@Override
	public double evaluateFaceIntegral(FT face, SystemShapeFunction<CT, FT, ET> shapeFunction1, SystemShapeFunction<CT, FT, ET> shapeFunction2)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
