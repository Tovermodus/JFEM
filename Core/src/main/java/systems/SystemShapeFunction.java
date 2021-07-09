package systems;

import basic.*;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Map;
import java.util.Set;

public class SystemShapeFunction<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>>
	extends SystemFunction
	implements ShapeFunction<CT,FT,ET,SystemShapeFunction<CT,FT,ET>,SystemValue, SystemGradient, SystemHessian>
{
	@Override
	public Set<CT> getCells()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Set<FT> getFaces()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public NodeFunctional<? extends Function<SystemValue, SystemGradient, SystemHessian>, SystemValue, SystemGradient, SystemHessian> getNodeFunctional()
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
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public SystemGradient gradientInCell(CoordinateVector pos, CT cell)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public SystemValue jumpInValue(FT face, CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public SystemGradient jumpInDerivative(FT face, CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public SystemValue averageInValue(FT face, CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public SystemGradient averageInDerivative(FT face, CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public SystemGradient normalAverageInValue(FT face, CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public SystemValue normalAverageInDerivative(FT face, CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Map<Integer, Double> prolongate(Set<SystemShapeFunction<CT, FT, ET>> refinedFunctions)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int compareTo(@NotNull SystemShapeFunction<CT, FT, ET> o)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
