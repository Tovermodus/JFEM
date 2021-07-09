package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.lang.reflect.InvocationTargetException;
import java.util.Set;

public class SingleComponentVectorShapeFunction<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>,
	ET extends Edge<CT,FT,ET>,
	ST extends FastEvaluatedScalarShapeFunction<CT,FT,ET,ST>, VST extends SingleComponentVectorShapeFunction<CT,
	FT,ET,ST,VST>> implements VectorShapeFunction<CT,FT,ET,
	VST>
{
	private final ST componentFunction;
	private final int component;
	private int globalIndex;
	
	public SingleComponentVectorShapeFunction(CT supportCell, int polynomialDegree, int localIndex,
	                             Class<ST> componentFunctionClass)
	{
		component = (int) (localIndex / Math.pow((polynomialDegree+1), supportCell.getDimension()));
		int componentLocalIndex = (int) (localIndex % Math.pow((polynomialDegree+1), supportCell.getDimension()));
		try
		{
			componentFunction =
				componentFunctionClass.getConstructor(supportCell.getClass(),
					Integer.TYPE,
					Integer.TYPE).newInstance(supportCell, polynomialDegree, componentLocalIndex);
		} catch (InstantiationException | IllegalAccessException | InvocationTargetException | NoSuchMethodException e)
		{
			throw new IllegalStateException();
		}
	}
	
	public ST getComponentFunction()
	{
		return componentFunction;
	}
	
	public int getComponent()
	{
		return component;
	}
	
	@Override
	public int getDomainDimension()
	{
		return componentFunction.getDomainDimension();
	}
	
	@Override
	public Set<CT> getCells()
	{
		return componentFunction.getCells();
	}
	
	@Override
	public Set<FT> getFaces()
	{
		return componentFunction.getFaces();
	}
	
	
	@Override
	public CoordinateVector valueInCell(CoordinateVector pos, CT cell)
	{
		return CoordinateVector.getUnitVector(getDomainDimension(), component).mul(componentFunction.fastValueInCell(pos,
			cell));
	}
	
	@Override
	public CoordinateMatrix gradientInCell(CoordinateVector pos, CT cell)
	{
		return componentFunction.gradientInCell(pos, cell).outer(CoordinateVector.getUnitVector(getDomainDimension(), component));
	}
	
	
	
	@Override
	public int compareTo(@NotNull VST o)
	{
		if(o.getComponent() < getComponent())
			return 1;
		else if(o.getComponent() > getComponent())
			return -1;
		else
			return componentFunction.compareTo(o.getComponentFunction());
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
	public double divergenceInCell(CoordinateVector pos, CT cell)
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
	public CoordinateVector curl(CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(getDomainDimension());
		CoordinateVector grad = componentFunction.gradient(pos);
		int compplus1 = (component+1)%3;
		int compplus2 = (component+2)%3;
		ret.set(grad.at(compplus2),compplus1);
		ret.set(-grad.at(compplus1),compplus2);
		return ret;
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
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	public void setGlobalIndex(int globalIndex)
	{
		this.globalIndex = globalIndex;
	}
}
