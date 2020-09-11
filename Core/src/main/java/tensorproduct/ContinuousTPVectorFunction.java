package tensorproduct;

import basic.SingleComponentVectorShapeFunction;

public class ContinuousTPVectorFunction extends SingleComponentVectorShapeFunction<TPCell,
	TPFace, ContinuousTPShapeFunction, ContinuousTPVectorFunction>
{
	public ContinuousTPVectorFunction(TPCell supportCell, int polynomialDegree, int localIndex,
	         Class<ContinuousTPShapeFunction> componentFunctionClass)
	{
		super(supportCell, polynomialDegree, localIndex, componentFunctionClass);
	}
}
