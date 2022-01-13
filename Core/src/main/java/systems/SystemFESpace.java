package systems;

import basic.Cell;
import basic.Face;
import basic.MatrixFESpace;
import com.google.common.collect.Multimap;
import linalg.CoordinateVector;
import linalg.MutableMatrix;
import linalg.MutableVector;

import java.util.Collection;
import java.util.List;
import java.util.Map;

public class SystemFESpace<CT extends Cell<CT, FT>, FT extends Face<CT, FT>
	>
	implements MatrixFESpace<CT, FT, SystemShapeFunction<CT, FT, ?>>
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
	public Map<Integer, SystemShapeFunction<CT, FT, ?>> getShapeFunctions()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public List<FT> getFaces()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Collection<SystemShapeFunction<CT, FT, ?>> getShapeFunctionsWithSupportOnCell(final CT cell)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Collection<SystemShapeFunction<CT, FT, ?>> getShapeFunctionsWithSupportOnFace(final FT face)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public List<CoordinateVector> generatePlotPoints(final int resolution)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Multimap<CT, SystemShapeFunction<CT, FT, ?>> getCellSupportMapping()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Multimap<FT, SystemShapeFunction<CT, FT, ?>> getFaceSupportMapping()
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
