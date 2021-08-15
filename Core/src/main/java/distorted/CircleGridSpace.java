package distorted;

import basic.AcceptsMatrixBoundaryValues;
import basic.Assembleable;
import basic.ShapeFunction;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Sets;
import distorted.geometry.CircleGrid;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.SparseMatrix;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public abstract class CircleGridSpace<ST extends ShapeFunction<DistortedCell, DistortedFace, valueT, gradientT, hessianT>,
	valueT, gradientT, hessianT>
	implements AcceptsMatrixBoundaryValues<DistortedCell, DistortedFace, ST, valueT, gradientT, hessianT>, Assembleable
{
	protected final CircleGrid grid;
	protected final HashMultimap<DistortedCell, ST> supportOnCell;
	protected final HashMultimap<DistortedFace, ST> supportOnFace;
	protected Set<ST> shapeFunctions;
	protected SparseMatrix systemMatrix;
	protected DenseVector rhs;
	Set<Integer> fixedNodes;
	
	@Override
	public Set<Integer> getFixedNodeIndices()
	{
		return fixedNodes;
	}
	
	public CircleGridSpace(final CoordinateVector center, final double radius, final int refinements)
	{
		supportOnCell = HashMultimap.create();
		supportOnFace = HashMultimap.create();
		grid = new CircleGrid(center, radius, refinements);
		fixedNodes = Sets.newConcurrentHashSet();
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
		systemMatrix = new SparseMatrix(shapeFunctions.size(), shapeFunctions.size());
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
	public List<DistortedCell> getCells()
	{
		return grid.cells.asList();
	}
	
	@Override
	public Map<Integer, ST> getShapeFunctions()
	{
		
		final Map<Integer, ST> functionNumbers = new TreeMap<>();
		for (final ST shapeFunction : shapeFunctions)
			functionNumbers.put(shapeFunction.getGlobalIndex(), shapeFunction);
		return functionNumbers;
	}
	
	@Override
	public List<DistortedFace> getFaces()
	{
		return grid.faces.asList();
	}
	
	@Override
	public Set<ST> getShapeFunctionsWithSupportOnCell(final DistortedCell cell)
	{
		return supportOnCell.get(cell);
	}
	
	@Override
	public Set<ST> getShapeFunctionsWithSupportOnFace(final DistortedFace face)
	{
		return supportOnFace.get(face);
	}
	
	@Override
	public List<CoordinateVector> generatePlotPoints(final int resolution)
	{
		return grid.generatePlotPoints(resolution);
	}
}
