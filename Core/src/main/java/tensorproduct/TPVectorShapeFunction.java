package tensorproduct;

import basic.*;
import linalg.*;
import org.jetbrains.annotations.NotNull;

import java.util.*;

public class TPVectorShapeFunction extends VectorShapeFunction<TPCell,TPFace,
	TPVectorShapeFunction> implements Comparable<TPVectorShapeFunction>
{
	TPShapeFunction componentFunction;
	int component;
	LagrangeVectorNodeFunctional nodeFunctional;
	public TPVectorShapeFunction(TPCell supportCell, int polynomialDegree, int localIndex)
	{
		component = (int) (localIndex / Math.pow((polynomialDegree+1), supportCell.getDimension()));
		int componentLocalIndex = (int) (localIndex % Math.pow((polynomialDegree+1), supportCell.getDimension()));
		componentFunction = new TPShapeFunction(supportCell,polynomialDegree, componentLocalIndex);
		CoordinateVector functional_point =
			((LagrangeNodeFunctional)(componentFunction.getNodeFunctional())).getPoint();
		nodeFunctional = new LagrangeVectorNodeFunctional(functional_point, component);
	}
	
	@Override
	public int getDomainDimension()
	{
		return componentFunction.getDomainDimension();
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
		throw new UnsupportedOperationException("TPShapefunctions do not live on faces");
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
	
	
	@Override
	public Map<Integer, Double> prolongate(Set<TPVectorShapeFunction> refinedFunctions)
	{
		Map<Integer, Double> ret = new HashMap<>();
		for(TPVectorShapeFunction shapeFunction:refinedFunctions)
		{
			ret.put(shapeFunction.getGlobalIndex(), shapeFunction.getNodeFunctional().evaluate(this));
		}
		return ret;
	}
	
	@Override
	public int compareTo(@NotNull TPVectorShapeFunction o)
	{
		if(o.component < component)
			return 1;
		else if(o.component > component)
			return -1;
		else
			return componentFunction.compareTo(o.componentFunction);
	}
	
	@Override
	public String toString()
	{
		return "Component: "+ component+ componentFunction.toString();
	}
	
	@Override
	public NodeFunctional<VectorFunction, CoordinateVector, CoordinateMatrix, Tensor> getNodeFunctional()
	{
		return nodeFunctional;
	}
	
	@Override
	public double divergenceInCell(CoordinateVector pos, TPCell cell)
	{
		return componentFunction.gradientInCell(pos, cell).at(component);
	}
	
	@Override
	public double divergence(CoordinateVector pos)
	{
		return divergenceInCell(pos, componentFunction.supportCell);
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
