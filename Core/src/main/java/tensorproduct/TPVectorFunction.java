package tensorproduct;

import basic.SingleComponentVectorShapeFunction;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TPVectorFunction extends SingleComponentVectorShapeFunction<TPCell,
	TPFace, TPShapeFunction>
{
	public TPVectorFunction(final TPCell supportCell, final int polynomialDegree, final int localIndex, final Class<TPShapeFunction> componentFunctionClass)
	{
		super(supportCell, polynomialDegree, localIndex, componentFunctionClass);
	}
	
	public TPVectorFunction(final TPCell supportCell, final int polynomialDegree, final int localIndex)
	{
		super(supportCell, polynomialDegree, localIndex, TPShapeFunction.class);
	}
}
