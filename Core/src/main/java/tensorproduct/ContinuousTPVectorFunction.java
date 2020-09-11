package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.SingleComponentVectorShapeFunction;
import basic.VectorNodeFunctional;
import linalg.CoordinateVector;

public class ContinuousTPVectorFunction extends SingleComponentVectorShapeFunction<TPCell,
	TPFace, ContinuousTPShapeFunction, ContinuousTPVectorFunction>
{
	public ContinuousTPVectorFunction(TPCell supportCell, int polynomialDegree, int localIndex,
	         Class<ContinuousTPShapeFunction> componentFunctionClass)
	{
		super(supportCell, polynomialDegree, localIndex, componentFunctionClass);
	}
	public CoordinateVector getNodeFunctionalPoint()
	{
		return ((LagrangeNodeFunctional)getNodeFunctional().getComponentNodeFunctional()).getPoint();
	}
}
