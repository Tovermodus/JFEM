package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.SingleComponentVectorShapeFunction;
import basic.VectorShapeFunctionWithReferenceShapeFunction;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class ContinuousTPVectorFunction extends SingleComponentVectorShapeFunction<TPCell,
	TPFace, ContinuousTPShapeFunction> implements VectorShapeFunctionWithReferenceShapeFunction<TPCell, TPFace>
{
	public ContinuousTPVectorFunction(TPCell supportCell, int polynomialDegree, int localIndex)
	{
		super(supportCell, polynomialDegree, localIndex, ContinuousTPShapeFunction.class);
	}
	public ContinuousTPVectorFunction(ContinuousTPShapeFunction function, int component)
	{
		super(function, component);
	}
	public CoordinateVector getNodeFunctionalPoint()
	{
		return ((LagrangeNodeFunctional)getNodeFunctional().getComponentNodeFunctional()).getPoint();
	}
	
	@Override
	public ContinuousTPVectorFunction createReferenceShapeFunctionRelativeTo(TPCell cell)
	{
		return new ContinuousTPVectorFunction(getComponentFunction().createReferenceShapeFunctionRelativeTo(cell)
			, getComponent());
	}
	
	@Override
	public ContinuousTPVectorFunction createReferenceShapeFunctionRelativeTo(TPFace face)
	{
		return new ContinuousTPVectorFunction(getComponentFunction().createReferenceShapeFunctionRelativeTo(face)
			, getComponent());
	}
	@Override
	public int hashCode()
	{
		return (17655+17*getComponentFunction().hashCode())*(93453+getComponent());
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(obj instanceof ContinuousTPVectorFunction)
		{
			ContinuousTPVectorFunction o = (ContinuousTPVectorFunction)obj;
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
