package tensorproduct;

import basic.SingleComponentVectorShapeFunction;

public class TPVectorFunction extends SingleComponentVectorShapeFunction<TPCell,
	TPFace,TPEdge, TPShapeFunction>
{
	public TPVectorFunction(TPCell supportCell, int polynomialDegree, int localIndex, Class<TPShapeFunction> componentFunctionClass)
	{
		super(supportCell, polynomialDegree, localIndex, componentFunctionClass);
	}
}
