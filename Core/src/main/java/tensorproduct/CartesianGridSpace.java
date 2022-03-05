package tensorproduct;

import basic.AcceptsMatrixBoundaryValues;
import basic.Assembleable;
import basic.ShapeFunction;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.IntCoordinates;
import linalg.SparseMatrix;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

public abstract class CartesianGridSpace<ST extends ShapeFunction<TPCell, TPFace, valueT, gradientT, hessianT>,
	valueT, gradientT, hessianT>
	implements AcceptsMatrixBoundaryValues<TPCell, TPFace, ST, valueT, gradientT, hessianT>, Assembleable
{
	public final CartesianGrid grid;
	private final HashMultimap<TPCell, ST> supportOnCell;
	private final HashMultimap<TPFace, ST> supportOnFace;
	protected Set<ST> shapeFunctions;
	protected SparseMatrix systemMatrix;
	protected DenseVector rhs;
	Set<Integer> fixedNodes;
	
	@Override
	public Set<Integer> getFixedNodeIndices()
	{
		return fixedNodes;
	}
	
	public CartesianGridSpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                          final List<Integer> cellsPerDimension)
	{
		this(startCoordinates, endCoordinates, new IntCoordinates(cellsPerDimension));
	}
	
	public CartesianGridSpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                          final IntCoordinates cellsPerDimension)
	{
		if (startCoordinates.getLength() != endCoordinates.getLength() ||
			startCoordinates.getLength() != cellsPerDimension.getDimension())
			throw new IllegalArgumentException();
		supportOnCell = HashMultimap.create();
		supportOnFace = HashMultimap.create();
		grid = new CartesianGrid(startCoordinates, endCoordinates, cellsPerDimension);
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
	public List<TPCell> getCells()
	{
		return grid.cells;
	}
	
	@Override
	public Map<Integer, ST> getShapeFunctionMap()
	{
		
		final Map<Integer, ST> functionNumbers = new Int2ObjectOpenHashMap<>();
		for (final ST shapeFunction : shapeFunctions)
			functionNumbers.put(shapeFunction.getGlobalIndex(), shapeFunction);
		return functionNumbers;
	}
	
	@Override
	public Collection<ST> getShapeFunctions()
	{
		return shapeFunctions;
	}
	
	@Override
	public List<TPFace> getFaces()
	{
		return grid.faces;
	}
	
	@Override
	public Set<ST> getShapeFunctionsWithSupportOnCell(final TPCell cell)
	{
		return supportOnCell.get(cell);
	}
	
	@Override
	public Set<ST> getShapeFunctionsWithSupportOnFace(final TPFace face)
	{
		return supportOnFace.get(face);
	}
	
	@Override
	public List<CoordinateVector> generatePlotPoints(final int resolution)
	{
		return grid.generatePlotPoints(resolution);
	}
	
	@Override
	public Multimap<TPCell, ST> getCellSupportMapping()
	{
		return supportOnCell;
	}
	
	@Override
	public Multimap<TPFace, ST> getFaceSupportMapping()
	{
		return supportOnFace;
	}
	
	public double getMaxDiam()
	{
		return grid.getMaxDiam();
	}
}
