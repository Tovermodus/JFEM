package systems;

import basic.*;
import com.google.common.collect.Sets;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class SystemShapeFunction<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>>
	extends SystemFunction
	implements ShapeFunction<CT,FT,ET,SystemShapeFunction<CT,FT,ET>,SystemValue, SystemGradient, SystemHessian>
{
	Set<CT> cells;
	Set<FT> faces;
	
	public SystemShapeFunction(ShapeFunction<CT,FT,ET, ?, ?, ?, ?>[] functions)
	{
		super(functions);
		cells = Arrays.stream(functions).map(ShapeFunction::getCells).reduce(Sets::union).orElse(new HashSet<>());
		faces = Arrays.stream(functions).map(ShapeFunction::getFaces).reduce(Sets::union).orElse(new HashSet<>());
	}
	
	@Override
	public Set<CT> getCells()
	{
		return cells;
	}
	
	@Override
	public Set<FT> getFaces()
	{
		return faces;
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
	@SuppressWarnings("unchecked")
	public SystemValue valueInCell(CoordinateVector pos, CT cell)
	{
		SystemValue ret = new SystemValue();
		for(int i = 0; i < functions.length; i++)
		{
			if(Double.class.isAssignableFrom(SystemParameters.getInstance().signatures[i].getValueT()))
				ret.setComponent(((ShapeFunction<CT,FT,ET,?,Double,CoordinateVector,
					CoordinateMatrix>) functions[i]).valueInCell(pos, cell), i);
			if(CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[i].getValueT()))
				ret.setComponent(((ShapeFunction<CT,FT,ET,?,CoordinateVector,
					CoordinateMatrix, CoordinateTensor>) functions[i]).valueInCell(pos, cell), i);
		}
		return ret;
	}
	
	@Override
	@SuppressWarnings("unchecked")
	public SystemGradient gradientInCell(CoordinateVector pos, CT cell)
	{
		SystemGradient ret = new SystemGradient(getDomainDimension());
		for(int i = 0; i < functions.length; i++)
		{
			if(Double.class.isAssignableFrom(SystemParameters.getInstance().signatures[i].getValueT()))
				ret.setComponent(((ShapeFunction<CT,FT,ET,?,Double,CoordinateVector,
					CoordinateMatrix>) functions[i]).gradientInCell(pos, cell), i);
			if(CoordinateVector.class.isAssignableFrom(SystemParameters.getInstance().signatures[i].getValueT()))
				ret.setComponent(((ShapeFunction<CT,FT,ET,?,CoordinateVector,
					CoordinateMatrix, CoordinateTensor>) functions[i]).gradientInCell(pos, cell), i);
		}
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
	public Map<Integer, Double> prolongate(Set<SystemShapeFunction<CT, FT, ET>> refinedFunctions)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	@SuppressWarnings("unchecked")
	public int compareTo(@NotNull SystemShapeFunction<CT, FT, ET> o)
	{
		int comp = 0;
		int i = 0;
		while (comp == 0)
		{
			comp = ((ShapeFunction<?, ?, ?, ShapeFunction, ?, ?, ?>) functions[i]).compareTo((ShapeFunction<?, ?, ?,
				ShapeFunction, ?, ?, ?>) (o.functions[i]));
			i++;
		}
		return comp;
	}
}
