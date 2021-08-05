package tensorproduct;

import basic.Assembleable;
import basic.MatrixFESpace;
import basic.ShapeFunction;
import com.google.common.collect.HashMultimap;
import linalg.*;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;

public abstract class CartesianGridSpace<ST extends ShapeFunction<TPCell, TPFace, ?, ?, ?>> implements MatrixFESpace<TPCell,
	TPFace, ST>, Assembleable
{
	protected final CartesianGrid grid;
	protected final HashMultimap<TPCell,ST> supportOnCell;
	protected final HashMultimap<TPFace,ST> supportOnFace;
	protected Set<ST> shapeFunctions;
	protected SparseMatrix systemMatrix;
	protected DenseVector rhs;
	private volatile int cellCounter = 0;
	private volatile int faceCounter = 0;
	
	public CartesianGridSpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                 List<Integer> cellsPerDimension)
	{
		if(startCoordinates.getLength() != endCoordinates.getLength()|| startCoordinates.getLength() != cellsPerDimension.size())
			throw new IllegalArgumentException();
		supportOnCell = HashMultimap.create();
		supportOnFace = HashMultimap.create();
		grid = new CartesianGrid(startCoordinates, endCoordinates, new IntCoordinates(cellsPerDimension));
	}
	@Override
	public void assembleCells()
	{
	
	}
	
	@Override
	public abstract void assembleFunctions(int polynomialDegree);
	
	@Override
	public void initializeSystemMatrix()
	{
		systemMatrix = new SparseMatrix(shapeFunctions.size(),shapeFunctions.size());
	}
	
	@Override
	public void initializeRhs()
	{
		rhs = new DenseVector(shapeFunctions.size());
	}
	
	@Override
	public DenseVector getRhs()
	{
		return rhs;
	}
	
	@Override
	public SparseMatrix getSystemMatrix()
	{
		return systemMatrix;
	}
	
	@Override
	public int getDimension()
	{
		return grid.dimension;
	}
	
	@Override
	public List<TPCell> getCells()
	{
		return grid.cells;
	}
	
	@Override
	public Map<Integer, ST> getShapeFunctions()
	{
		
		Map<Integer, ST> functionNumbers = new TreeMap<>();
		for(ST shapeFunction:shapeFunctions)
			functionNumbers.put(shapeFunction.getGlobalIndex(), shapeFunction);
		return functionNumbers;
	}
	
	@Override
	public List<TPFace> getFaces()
	{
		return grid.faces;
	}
	
	@Override
	public Set<ST> getShapeFunctionsWithSupportOnCell(TPCell cell)
	{
		return supportOnCell.get(cell);
	}
	
	@Override
	public Set<ST> getShapeFunctionsWithSupportOnFace(TPFace face)
	{
		return supportOnFace.get(face);
	}
	
	@Override
	public List<CoordinateVector> generatePlotPoints(int resolution)
	{
		return grid.generatePlotPoints(resolution);
	}
	
	
	@Override
	public synchronized int increaseCellCounter()
	{
		return cellCounter++;
	}
	
	@Override
	public synchronized int increaseFaceCounter()
	{
		return faceCounter++;
	}
	
}
