package distorted;

import basic.SingleComponentVectorShapeFunction;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;

public class DistortedVectorShapeFunction extends SingleComponentVectorShapeFunction<DistortedCell, DistortedFace,
	DistortedShapeFunction, DistortedVectorShapeFunction>

{
	public DistortedVectorShapeFunction(final DistortedCell supportCell, final int polynomialDegree, final int localIndex)
	{
		super(supportCell, polynomialDegree, localIndex, DistortedShapeFunction.class);
	}
	
	public DistortedVectorShapeFunction(final DistortedShapeFunction componentFunction, final int component)
	{
		super(componentFunction, component);
	}
}
