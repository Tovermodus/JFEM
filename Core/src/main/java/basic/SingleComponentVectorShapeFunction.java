package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

import java.lang.reflect.InvocationTargetException;
import java.util.Objects;
import java.util.Set;

public class SingleComponentVectorShapeFunction<CT extends Cell<CT, FT>, FT extends Face<CT, FT>, CST extends FastEvaluatedScalarShapeFunction<CT, FT> & Comparable<CST>, VT extends SingleComponentVectorShapeFunction<CT, FT, CST, VT>>
	implements VectorShapeFunction<CT, FT>, Comparable<VT>
{
	private final CST componentFunction;
	private final int component;
	private int globalIndex = -1;
	private final int dimension;
	
	public SingleComponentVectorShapeFunction(final CT supportCell,
	                                          final int polynomialDegree,
	                                          final int localIndex,
	                                          final Class<CST> componentFunctionClass)
	{
		component = (int) (localIndex / Math.pow((polynomialDegree + 1), supportCell.getDimension()));
		dimension = supportCell.getDimension();
		final int componentLocalIndex = (int) (localIndex % Math.pow((polynomialDegree + 1),
		                                                             supportCell.getDimension()));
		try
		{
			componentFunction = componentFunctionClass
				.getConstructor(supportCell.getClass(), Integer.TYPE, Integer.TYPE)
				.newInstance(supportCell, polynomialDegree, componentLocalIndex);
		} catch (final InstantiationException | IllegalAccessException | InvocationTargetException | NoSuchMethodException e)
		{
			e.printStackTrace();
			throw new IllegalStateException();
		}
	}
	
	public SingleComponentVectorShapeFunction(final CST componentFunction, final int component)
	{
		dimension = componentFunction.getDomainDimension();
		this.component = component;
		this.componentFunction = componentFunction;
	}
	
	public CST getComponentFunction()
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
		return dimension;
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
	public CoordinateVector valueInCell(final CoordinateVector pos, final CT cell)
	{
		return CoordinateVector
			.getUnitVector(getDomainDimension(), component)
			.mul(componentFunction.fastValueInCell(pos, cell));
	}
	
	@Override
	public CoordinateMatrix gradientInCell(final CoordinateVector pos, final CT cell)
	{
		return componentFunction
			.gradientInCell(pos, cell)
			.outer(CoordinateVector.getUnitVector(getDomainDimension(), component));
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
	
	@Override
	public double divergenceInCell(final CoordinateVector pos, final CT cell)
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
	public CoordinateVector curl(final CoordinateVector pos)
	{
		final CoordinateVector ret = new CoordinateVector(getDomainDimension());
		final CoordinateVector grad = componentFunction.gradient(pos);
		final int compplus1 = (component + 1) % 3;
		final int compplus2 = (component + 2) % 3;
		ret.set(grad.at(compplus2), compplus1);
		ret.set(-grad.at(compplus1), compplus2);
		return ret;
	}
	
	@Override
	public CoordinateVector value(final CoordinateVector pos)
	{
		return CoordinateVector
			.getUnitVector(getDomainDimension(), component)
			.mul(componentFunction.fastValue(pos));
	}
	
	@Override
	public CoordinateMatrix gradient(final CoordinateVector pos)
	{
		return componentFunction.gradient(pos)
		                        .outer(CoordinateVector.getUnitVector(getDomainDimension(), component));
	}
	
	@Override
	public CoordinateMatrix jacobian(final CoordinateVector pos)
	{
		return CoordinateVector.getUnitVector(getDomainDimension(), component)
		                       .outer(componentFunction.gradient(pos));
	}
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	public void setGlobalIndex(final int globalIndex)
	{
		this.globalIndex = globalIndex;
	}
	
	@Override
	public int compareTo(final VT o)
	{
		if (o.getComponent() < getComponent()) return 1;
		else if (o.getComponent() > getComponent()) return -1;
		else return (componentFunction).compareTo(o.getComponentFunction());
	}
	
	@Override
	public boolean equals(final Object o)
	{
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		final SingleComponentVectorShapeFunction<?, ?, ?, ?> that
			= (SingleComponentVectorShapeFunction<?, ?, ?, ?>) o;
		return component == that.component && dimension == that.dimension && Objects
			.equals(componentFunction, that.componentFunction);
	}
	
	@Override
	public int hashCode()
	{
		return Objects.hash(componentFunction, component, dimension);
	}
}
