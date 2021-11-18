package systems;

import basic.Cell;
import basic.Face;
import basic.ShapeFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

public class SystemShapeFunction<CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
	CST extends ShapeFunction<CT, FT, ?, ?, ?>>
	extends SystemFunction
	implements ShapeFunction<CT, FT, SystemValue, SystemGradient,
	SystemHessian>, Comparable<SystemShapeFunction<CT, FT, CST>>
{
	private final CST function;
	
	public SystemShapeFunction(final CST function, final int component)
	{
		super(function, component);
		this.function = function;
	}
	
	@Override
	public Collection<CT> getCells()
	{
		return function.getCells();
	}
	
	@Override
	public Collection<FT> getFaces()
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
	public SystemValue valueInCell(final CoordinateVector pos, final CT cell)
	{
		final SystemValue ret = new SystemValue();
		if (Double.class.isAssignableFrom(SystemParameters.getInstance().signatures[mainComponent].getValueT()))
			ret.setComponent((Double) function.valueInCell(pos, cell), mainComponent);
		if (CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[mainComponent].getValueT()))
			ret.setComponent((CoordinateVector) function.valueInCell(pos, cell), mainComponent);
		return ret;
	}
	
	@Override
	public SystemGradient gradientInCell(final CoordinateVector pos, final CT cell)
	{
		final SystemGradient ret = new SystemGradient(getDomainDimension());
		if (CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[mainComponent].getGradientT()))
			ret.setComponent((CoordinateVector) function.gradientInCell(pos, cell), mainComponent);
		if (CoordinateMatrix.class.isAssignableFrom(SystemParameters.getInstance().signatures[mainComponent].getGradientT()))
			ret.setComponent((CoordinateMatrix) function.gradientInCell(pos, cell), mainComponent);
		return ret;
	}
	
	@Override
	public SystemValue jumpInValue(final FT face, final CoordinateVector pos)
	{
		return new SystemValue(valueInCell(pos, face.getNormalUpstreamCell()).sub(valueInCell(pos,
		                                                                                      face.getNormalDownstreamCell())),
		                       true);
	}
	
	@Override
	public SystemGradient jumpInDerivative(final FT face, final CoordinateVector pos)
	{
		return new SystemGradient(gradientInCell(pos, face.getNormalUpstreamCell()).sub(gradientInCell(pos,
		                                                                                               face.getNormalDownstreamCell())),
		                          true);
	}
	
	@Override
	public SystemValue averageInValue(final FT face, final CoordinateVector pos)
	{
		return new SystemValue(valueInCell(pos, face.getNormalUpstreamCell()).add(valueInCell(pos,
		                                                                                      face.getNormalDownstreamCell()))
		                                                                     .mul(0.5), true);
	}
	
	@Override
	public SystemGradient averageInDerivative(final FT face, final CoordinateVector pos)
	{
		return new SystemGradient(gradientInCell(pos, face.getNormalUpstreamCell()).add(
			gradientInCell(pos, face.getNormalDownstreamCell()))
		                                                                           .mul(0.5), true);
	}
	
	@Override
	public SystemGradient normalAverageInValue(final FT face, final CoordinateVector pos)
	{
		return new SystemGradient(face.getNormal()
		                              .value(pos)
		                              .outer(jumpInValue(face, pos).mul(0.5)), true);
	}
	
	@Override
	public SystemValue normalAverageInDerivative(final FT face, final CoordinateVector pos)
	{
		return new SystemValue(jumpInDerivative(face, pos).mvMul(face.getNormal()
		                                                             .value(pos)), true);
	}
	
	@Override
	public <ST extends ShapeFunction<CT, FT, SystemValue, SystemGradient, SystemHessian>> Map<Integer, Double> prolongate(
		final Set<ST> refinedFunctions)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public CST getComponentFunction(final int component)
	{
		return function;
	}
	
	@Override
	public String toString()
	{
		return "component: " + mainComponent + " " + function.getClass() + " function: " + function;
	}
	
	@Override
	@SuppressWarnings("unchecked")
	public int compareTo(final SystemShapeFunction<CT, FT, CST> o)
	{
		if (o.mainComponent == mainComponent)
			return ((Comparable<CST>) function).compareTo(o.function);
		return Integer.compare(mainComponent, o.mainComponent);
	}
}
