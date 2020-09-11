package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Matrix;
import linalg.Tensor;
import org.jetbrains.annotations.NotNull;
import tensorproduct.TPCell;
import tensorproduct.TPFace;
import tensorproduct.TPShapeFunction;

import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.Set;

public class SingleComponentVectorShapeFunction<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	ST extends ScalarShapeFunction<CT,FT,ST>, VST extends SingleComponentVectorShapeFunction<CT,FT,ST,VST>> extends VectorShapeFunction<CT,FT,
	VST>
{
	private final ST componentFunction;
	private final int component;
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
	public void addFace(FT face)
	{
		componentFunction.addFace(face);
	}
	
	@Override
	public void addCell(CT cell)
	{
		componentFunction.addCell(cell);
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
