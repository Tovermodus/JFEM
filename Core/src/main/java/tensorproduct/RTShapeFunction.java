package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.VectorNodeFunctional;
import basic.VectorShapeFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.Set;

public class RTShapeFunction
	implements VectorShapeFunction<TPCell, TPFace>, Comparable<RTShapeFunction>
{
	private final RTComponentFunction componentFunction;
	private final int component;
	
	public RTShapeFunction(final TPCell supportCell, final int polynomialDegree, final int localIndex)
	{
		component = (int) (localIndex / ((polynomialDegree + 2) * Math.pow((polynomialDegree + 1),
		                                                                   supportCell.getDimension() - 1)));
		final int componentLocalIndex = (int) (localIndex % ((polynomialDegree + 2) * Math.pow((polynomialDegree + 1),
		                                                                                       supportCell.getDimension() - 1)));
		componentFunction = new RTComponentFunction(supportCell, polynomialDegree, componentLocalIndex,
		                                            component);
	}
	
	@Override
	public Set<TPCell> getCells()
	{
		return componentFunction.getCells();
	}
	
	@Override
	public Set<TPFace> getFaces()
	{
		return componentFunction.getFaces();
	}
	
	@Override
	public CoordinateVector valueInCell(final CoordinateVector pos, final TPCell cell)
	{
		return CoordinateVector.getUnitVector(getDomainDimension(), component)
		                       .mul(componentFunction.fastValueInCell(pos,
		                                                              cell));
	}
	
	@Override
	public CoordinateMatrix gradientInCell(final CoordinateVector pos, final TPCell cell)
	{
		return componentFunction.gradientInCell(pos, cell)
		                        .outer(CoordinateVector.getUnitVector(getDomainDimension(), component));
	}
	
	public CoordinateVector getNodeFunctionalPoint()
	{
		return ((LagrangeNodeFunctional) getNodeFunctional().getComponentNodeFunctional()).getPoint();
	}
	
	@Override
	public int compareTo(@NotNull final RTShapeFunction o)
	{
		if (o.getComponent() < getComponent())
			return 1;
		else if (o.getComponent() > getComponent())
			return -1;
		else
			return componentFunction.compareTo(o.getComponentFunction());
	}
	
	@Override
	public int getDomainDimension()
	{
		return componentFunction.getDomainDimension();
	}
	
	public int getComponent()
	{
		return component;
	}
	
	public RTComponentFunction getComponentFunction()
	{
		return componentFunction;
	}
	
	@Override
	public String toString()
	{
		return "Component: " + component + componentFunction.toString();
	}
	
	@Override
	public VectorNodeFunctional getNodeFunctional()
	{
		return new VectorNodeFunctional(component, componentFunction.getNodeFunctional());
	}
	
	public void setGlobalIndex(final int index)
	{
		componentFunction.setGlobalIndex(index);
	}
	
	@Override
	public int getGlobalIndex()
	{
		return componentFunction.getGlobalIndex();
	}
	
	@Override
	public double divergenceInCell(final CoordinateVector pos, final TPCell cell)
	{
		return componentFunction.gradientInCell(pos, cell)
		                        .at(component);
	}
	
	@Override
	public int getRangeDimension()
	{
		return getDomainDimension();
	}
	
	@Override
	public double divergence(final CoordinateVector pos)
	{
		return componentFunction.gradient(pos)
		                        .at(component);
	}
	
	@Override
	public CoordinateVector value(final CoordinateVector pos)
	{
		return CoordinateVector.getUnitVector(getDomainDimension(), component)
		                       .mul(componentFunction.value(pos));
	}
	
	@Override
	public CoordinateMatrix gradient(final CoordinateVector pos)
	{
		return componentFunction.gradient(pos)
		                        .outer(CoordinateVector.getUnitVector(getDomainDimension(),
		                                                              component));
	}
}
