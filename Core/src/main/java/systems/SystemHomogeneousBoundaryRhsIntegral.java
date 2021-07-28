package systems;

import basic.BoundaryRightHandSideIntegral;
import basic.Face;
import basic.ShapeFunction;

public class SystemHomogeneousBoundaryRhsIntegral<FT extends Face<?,FT,?>, ST extends ShapeFunction<?,FT,?,?,?,?>>
	extends BoundaryRightHandSideIntegral<FT, SystemShapeFunction<?,FT,?,ST>>
{
	private final BoundaryRightHandSideIntegral<FT, ST> faceIntegral;
	private final int component;
	
	public SystemHomogeneousBoundaryRhsIntegral(BoundaryRightHandSideIntegral<FT, ST> faceIntegral, int component)
	{
		this.faceIntegral = faceIntegral;
		this.component = component;
	}
	
	
	@Override
	public double evaluateBoundaryRightHandSideIntegral(FT face, SystemShapeFunction<?, FT, ?, ST> shapeFunction)
	{
		if(shapeFunction.mainComponent == component)
		{
			return this.faceIntegral.evaluateBoundaryRightHandSideIntegral(face,
				shapeFunction.getComponentFunction(component));
		}
		return 0;
	}
	
}
