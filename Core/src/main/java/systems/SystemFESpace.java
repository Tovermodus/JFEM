package systems;

import basic.Cell;
import basic.Edge;
import basic.Face;
import basic.MatrixFESpace;
import linalg.CoordinateVector;
import linalg.MutableMatrix;
import linalg.MutableVector;
import tensorproduct.TPShapeFunction;

import java.util.Collection;
import java.util.List;
import java.util.Map;

public class SystemFESpace<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>,
	ET extends Edge<CT,FT,ET>> implements MatrixFESpace<CT, FT, ET, SystemShapeFunction<CT, FT, ET, ?>, SystemValue, SystemGradient, SystemHessian>
{
	@Override
	public int getDimension()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public List<CT> getCells()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Map<Integer, SystemShapeFunction<CT, FT, ET, ?>> getShapeFunctions()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public List<FT> getFaces()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Collection<SystemShapeFunction<CT, FT, ET, ?>> getShapeFunctionsWithSupportOnCell(CT cell)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Collection<SystemShapeFunction<CT, FT, ET, ?>> getShapeFunctionsWithSupportOnFace(FT face)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public List<CoordinateVector> generatePlotPoints(int resolution)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public void initializeSystemMatrix()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public void initializeRhs()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public MutableVector getRhs()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public MutableMatrix getSystemMatrix()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
