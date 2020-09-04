package tensorproduct;

import basic.SingleComponentVectorShapeFunction;

public class TPVectorFunction extends SingleComponentVectorShapeFunction<TPCell,
	TPFace, TPShapeFunction, TPVectorFunction>
{
	public TPVectorFunction(TPCell supportCell, int polynomialDegree, int localIndex, Class<TPShapeFunction> componentFunctionClass)
	{
		super(supportCell, polynomialDegree, localIndex, componentFunctionClass);
	}
}
