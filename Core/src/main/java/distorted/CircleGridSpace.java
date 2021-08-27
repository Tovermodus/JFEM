package distorted;

import basic.AcceptsMatrixBoundaryValues;
import basic.Assembleable;
import basic.ScalarFunction;
import basic.ShapeFunction;
import com.google.common.collect.Sets;
import distorted.geometry.CircleGrid;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.SparseMatrix;

import java.util.*;

public abstract class CircleGridSpace<ST extends ShapeFunction<DistortedCell, DistortedFace, valueT, gradientT, hessianT> & Comparable<ST>, valueT, gradientT, hessianT> implements AcceptsMatrixBoundaryValues<DistortedCell, DistortedFace, ST, valueT, gradientT, hessianT>, Assembleable
{
	protected final CircleGrid grid;
	private final HashMap<DistortedCell, TreeSet<ST>> supportOnCell;
	private final HashMap<DistortedFace, TreeSet<ST>> supportOnFace;
	protected TreeSet<ST> shapeFunctions;
	protected SparseMatrix systemMatrix;
	protected DenseVector rhs;
	Set<Integer> fixedNodes;
	public final double radius;
	public final CoordinateVector center;
	
	@Override
	public Set<Integer> getFixedNodeIndices()
	{
		return fixedNodes;
	}
	
	public CircleGridSpace(final CoordinateVector center, final double radius, final int refinements)
	{
		supportOnCell = new HashMap<>();
		supportOnFace = new HashMap<>();
		this.radius = radius;
		this.center = center;
		grid = new CircleGrid(center, radius, refinements);
		fixedNodes = Sets.newConcurrentHashSet();
	}
	
	protected void addFunctionToCell(final ST function, final DistortedCell cell)
	{
		if (!supportOnCell.containsKey(cell)) supportOnCell.put(cell, new TreeSet<>());
		supportOnCell.get(cell).add(function);
	}
	
	protected void addFunctionToFace(final ST function, final DistortedFace face)
	{
		if (!supportOnFace.containsKey(face)) supportOnFace.put(face, new TreeSet<>());
		supportOnFace.get(face).add(function);
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
	
	public ScalarFunction getIndicatorFunction()
	{
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return getDimension();
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				final boolean isInSide = getCells().stream().anyMatch(c -> c.isInCellPrecise(pos));
				if (isInSide) return 1.;
				else return 0.;
			}
		};
	}
}
