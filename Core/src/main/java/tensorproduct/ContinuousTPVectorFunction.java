package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.SingleComponentVectorShapeFunction;
import basic.VectorShapeFunctionWithReferenceShapeFunction;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class ContinuousTPVectorFunction extends SingleComponentVectorShapeFunction<TPCell,
	TPFace, ContinuousTPShapeFunction, ContinuousTPVectorFunction> implements VectorShapeFunctionWithReferenceShapeFunction<TPCell, TPFace>
{
	public ContinuousTPVectorFunction(final TPCell supportCell, final int polynomialDegree, final int localIndex)
	{
		super(supportCell, polynomialDegree, localIndex, ContinuousTPShapeFunction.class);
	}
	
	public ContinuousTPVectorFunction(final ContinuousTPShapeFunction function, final int component)
	{
		super(function, component);
	}
	
	public CoordinateVector getNodeFunctionalPoint()
	{
		return ((LagrangeNodeFunctional) getNodeFunctional().getComponentNodeFunctional()).getPoint();
	}
	
	@Override
	public ContinuousTPVectorFunction createReferenceShapeFunctionRelativeTo(final TPCell cell)
	{
		return new ContinuousTPVectorFunction(
			getComponentFunction().createReferenceShapeFunctionRelativeTo(cell)
			, getComponent());
	}
	
	@Override
	public ContinuousTPVectorFunction createReferenceShapeFunctionRelativeTo(final TPFace face)
	{
		return new ContinuousTPVectorFunction(
			getComponentFunction().createReferenceShapeFunctionRelativeTo(face)
			, getComponent());
	}
	
	@Override
	public int hashCode()
	{
		return (17655 + 17 * getComponentFunction().hashCode()) * (93453 + getComponent());
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (obj instanceof ContinuousTPVectorFunction)
		{
			final ContinuousTPVectorFunction o = (ContinuousTPVectorFunction) obj;
			if (o.getComponent() < getComponent())
				return false;
			else if (o.getComponent() > getComponent())
				return false;
			else
				return getComponentFunction().compareTo(o.getComponentFunction()) == 0;
		}
		return false;
	}
}
