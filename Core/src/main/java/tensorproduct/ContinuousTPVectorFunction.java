package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.SingleComponentVectorShapeFunction;
import basic.VectorShapeFunctionWithReferenceShapeFunction;
import linalg.CoordinateVector;

public class ContinuousTPVectorFunction extends SingleComponentVectorShapeFunction<TPCell,
	TPFace, TPEdge, ContinuousTPShapeFunction> implements VectorShapeFunctionWithReferenceShapeFunction<TPCell, TPFace, TPEdge>
{
	public ContinuousTPVectorFunction(TPCell supportCell, int polynomialDegree, int localIndex,
	                                  Class<ContinuousTPShapeFunction> componentFunctionClass)
	{
		super(supportCell, polynomialDegree, localIndex, componentFunctionClass);
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
