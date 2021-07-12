package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.ShapeFunction;
import basic.VectorNodeFunctional;
import basic.VectorShapeFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.Set;

public class NodalNedelecShapeFunction implements VectorShapeFunction<TPCell, TPFace, TPEdge>
{
	private final NodalNedelecComponentFunction componentFunction;
	private final int component;
	public NodalNedelecShapeFunction(TPCell supportCell, int polynomialDegree, int localIndex)
	{
		component = (int) (localIndex / ((polynomialDegree+1)*Math.pow((polynomialDegree+2),
			supportCell.getDimension()-1)));
		int componentLocalIndex = (int) (localIndex % ((polynomialDegree+1)*Math.pow((polynomialDegree+2),
			supportCell.getDimension()-1)));
		componentFunction = new NodalNedelecComponentFunction(supportCell, polynomialDegree, componentLocalIndex,
				component);
	}
	public void setGlobalIndex(int index)
	{
		componentFunction.setGlobalIndex(index);
	}
	@Override
	public int getGlobalIndex()
	{
		return componentFunction.getGlobalIndex();
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
	public CoordinateVector curl(CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(getDomainDimension());
		int compplus1 = (component + 1) % 3;
		int compplus2 = (component + 2) % 3;
		List<? extends Function1D> function1Ds = null;
		for (TPCell cell : componentFunction.getCells())
			if (cell.isInCell(pos))
			{
				function1Ds = componentFunction.get1DFunctionsInCell(cell);
				break;
			}
		if (function1Ds == null)
			return ret;
		double component1 = 1;
		double component2 = 1;
		for (int j = 0; j < pos.getLength(); j++)
		{
			if (compplus1 == j)
				component1 *= function1Ds.get(j).derivative(pos.at(j));
			else
				component1 *= function1Ds.get(j).value(pos.at(j));
			if (compplus2 == j)
				component2 *= function1Ds.get(j).derivative(pos.at(j));
			else
				component2 *= function1Ds.get(j).value(pos.at(j));
		}
		ret.set(component2, compplus1);
		ret.set(-component1, compplus2);
		return ret;
	}
	public int compareTo(@NotNull NodalNedelecShapeFunction o)
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
	
	public NodalNedelecComponentFunction getComponentFunction()
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
	public int getRangeDimension()
	{
		return getDomainDimension();
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
