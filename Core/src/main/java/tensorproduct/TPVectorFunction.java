package tensorproduct;

import basic.SingleComponentVectorShapeFunction;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TPVectorFunction extends SingleComponentVectorShapeFunction<TPCell,
	TPFace, TPShapeFunction>
{
	public TPVectorFunction(TPCell supportCell, int polynomialDegree, int localIndex, Class<TPShapeFunction> componentFunctionClass)
	{
		super(supportCell, polynomialDegree, localIndex, componentFunctionClass);
	}
}
