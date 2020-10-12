package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.VectorNodeFunctional;
import basic.VectorShapeFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Set;

public class RTShapeFunction extends VectorShapeFunction<TPCell, TPFace,TPEdge, RTShapeFunction>
{
	private final RTComponentFunction componentFunction;
	private final int component;
	public RTShapeFunction(TPCell supportCell, int polynomialDegree, int localIndex)
	{
		component = (int) (localIndex / ((polynomialDegree+2)*Math.pow((polynomialDegree+1),
			supportCell.getDimension()-1)));
		int componentLocalIndex = (int) (localIndex % ((polynomialDegree+2)*Math.pow((polynomialDegree+1),
			supportCell.getDimension()-1)));
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
	public void addFace(TPFace face)
	{
		componentFunction.addFace(face);
	}
	
	@Override
	public void addCell(TPCell cell)
	{
		componentFunction.addCell(cell);
	}
	@Override
	public CoordinateVector valueInCell(CoordinateVector pos, TPCell cell)
	{
		return CoordinateVector.getUnitVector(getDomainDimension(), component).mul(componentFunction.fastValueInCell(pos,
			cell));
	}
	
	@Override
	public CoordinateMatrix gradientInCell(CoordinateVector pos, TPCell cell)
	{
		return componentFunction.gradientInCell(pos, cell).outer(CoordinateVector.getUnitVector(getDomainDimension(), component));
	}
	
	
	public CoordinateVector getNodeFunctionalPoint()
	{
		return ((LagrangeNodeFunctional)getNodeFunctional().getComponentNodeFunctional()).getPoint();
	}
	
	@Override
	public int compareTo(@NotNull RTShapeFunction o)
	{
		if(o.getComponent() < getComponent())
			return 1;
		else if(o.getComponent() > getComponent())
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
		return "Component: "+ component+ componentFunction.toString();
	}
	
	@Override
	public VectorNodeFunctional getNodeFunctional()
	{
		return new VectorNodeFunctional(component, componentFunction.getNodeFunctional());
	}
	
	@Override
	public double divergenceInCell(CoordinateVector pos, TPCell cell)
	{
		return componentFunction.gradientInCell(pos, cell).at(component);
	}
	
	@Override
	public double divergence(CoordinateVector pos)
	{
		return componentFunction.gradient(pos).at(component);
	}
	
	@Override
	public CoordinateVector value(CoordinateVector pos)
	{
		return CoordinateVector.getUnitVector(getDomainDimension(), component).mul(componentFunction.fastValue(pos));
	}
	
	@Override
	public CoordinateMatrix gradient(CoordinateVector pos)
	{
		return componentFunction.gradient(pos).outer(CoordinateVector.getUnitVector(getDomainDimension(),
			component));
	}
}
