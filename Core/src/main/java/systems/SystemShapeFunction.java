package systems;

import basic.*;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Map;
import java.util.Set;

public class SystemShapeFunction<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>,
	CST extends ShapeFunction<CT,FT,ET,?, ?, ?>>
	extends SystemFunction
	implements ShapeFunction<CT, FT, ET, SystemValue, SystemGradient,
	SystemHessian>
{
	private final CST function;
	public SystemShapeFunction(CST function, int component)
	{
		super(function, component);
		this.function = function;
	}
	@Override
	public Set<CT> getCells()
	{
		return function.getCells();
	}
	
	@Override
	public Set<FT> getFaces()
	{
		return function.getFaces();
	}
	
	@Override
	public SystemNodeFunctional getNodeFunctional()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int getGlobalIndex()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public SystemValue valueInCell(CoordinateVector pos, CT cell)
	{
		SystemValue ret = new SystemValue();
		if(Double.class.isAssignableFrom(SystemParameters.getInstance().signatures[mainComponent].getValueT()))
			ret.setComponent((Double)function.valueInCell(pos, cell), mainComponent);
		if(CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[mainComponent].getValueT()))
			ret.setComponent((CoordinateVector) function.valueInCell(pos, cell), mainComponent);
		return ret;
	}
	
	@Override
	public SystemGradient gradientInCell(CoordinateVector pos, CT cell)
	{
		SystemGradient ret = new SystemGradient(getDomainDimension());
		if(CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[mainComponent].getGradientT()))
			ret.setComponent((CoordinateVector) function.gradientInCell(pos, cell), mainComponent);
		if(CoordinateMatrix.class.isAssignableFrom(SystemParameters.getInstance().signatures[mainComponent].getGradientT()))
			ret.setComponent((CoordinateMatrix) function.gradientInCell(pos, cell), mainComponent);
		return ret;
	}
	
	@Override
	public SystemValue jumpInValue(FT face, CoordinateVector pos)
	{
		return new SystemValue(valueInCell(pos, face.getNormalUpstreamCell()).sub(valueInCell(pos,
			face.getNormalDownstreamCell())), true);
	}
	
	@Override
	public SystemGradient jumpInDerivative(FT face, CoordinateVector pos)
	{
		return new SystemGradient(gradientInCell(pos, face.getNormalUpstreamCell()).sub(gradientInCell(pos,
			face.getNormalDownstreamCell())), true);
	}
	
	@Override
	public SystemValue averageInValue(FT face, CoordinateVector pos)
	{
		return new SystemValue(valueInCell(pos,face.getNormalUpstreamCell()).add(valueInCell(pos,
			face.getNormalDownstreamCell())).mul(0.5), true);
	}
	
	@Override
	public SystemGradient averageInDerivative(FT face, CoordinateVector pos)
	{
		return new SystemGradient(gradientInCell(pos,face.getNormalUpstreamCell()).add(
			gradientInCell(pos,face.getNormalDownstreamCell())).mul(0.5), true);
	}
	
	@Override
	public SystemGradient normalAverageInValue(FT face, CoordinateVector pos)
	{
		return new SystemGradient(face.getNormal().value(pos).outer(jumpInValue(face,pos).mul(0.5)), true);
	}
	
	@Override
	public SystemValue normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		return new SystemValue(jumpInDerivative(face, pos).mvMul(face.getNormal().value(pos)), true);
	}
	
	@Override
	public <ST extends ShapeFunction<CT, FT, ET, SystemValue, SystemGradient, SystemHessian>> Map<Integer, Double> prolongate(Set<ST> refinedFunctions)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	
	@Override
	public CST getComponentFunction(int component)
	{
		return function;
	}
	
	@Override
	public String toString()
	{
		return "component: " + mainComponent + " " + function.getClass() + " function: " + function;
	}
	
	public int compareTo(SystemShapeFunction<CT,FT,ET,CST> o)
	{
		if(o.mainComponent == mainComponent)
			return function.compareTo(o.function);
		return Integer.compare(mainComponent, o.mainComponent);
	}
	
}
