package systems;

import basic.*;
import tensorproduct.QuadratureRule1D;

public class SystemHomogeneousFaceIntegral<FT extends Face<?,FT>, ST extends ShapeFunction<?,FT, ?,?,?>>
	extends FaceIntegral<FT, SystemShapeFunction<?,FT, ST>>
{
	private final FaceIntegral<FT, ST> faceIntegral;
	private final int component;
	
	public SystemHomogeneousFaceIntegral(FaceIntegral<FT, ST> faceIntegral, int component, QuadratureRule1D quadratureRule1D)
	{
		super("homogeneous", quadratureRule1D);
		this.faceIntegral = faceIntegral;
		this.component = component;
	}
	public SystemHomogeneousFaceIntegral(FaceIntegral<FT, ST> faceIntegral, int component)
	{
		super("homogeneous", QuadratureRule1D.Gauss5);
		this.faceIntegral = faceIntegral;
		this.component = component;
	}
	
	
	@Override
	public double evaluateFaceIntegral(FT face, SystemShapeFunction<?, FT, ST> shapeFunction1,
	                                   SystemShapeFunction<?, FT, ST> shapeFunction2)
	{
		if(shapeFunction1.mainComponent == component && shapeFunction2.mainComponent == component)
		{
			return this.faceIntegral.evaluateFaceIntegral(face,
				shapeFunction1.getComponentFunction(component),
				shapeFunction2.getComponentFunction(component));
		}
		return 0;
	}
	
}
