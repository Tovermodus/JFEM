package distorted;

import basic.AcceptsMatrixBoundaryValues;
import basic.Assembleable;
import basic.ScalarFunction;
import basic.ShapeFunction;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.google.common.collect.TreeMultimap;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import distorted.geometry.DistortedGrid;
import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.SparseMatrix;

import java.util.*;

public abstract class DistortedGridSpace<ST extends ShapeFunction<DistortedCell,
	DistortedFace, valueT,
	gradientT,
	hessianT> & Comparable<ST>, valueT, gradientT, hessianT>
	implements AcceptsMatrixBoundaryValues<DistortedCell, DistortedFace, ST, valueT, gradientT, hessianT>,
	Assembleable
{
	protected final DistortedGrid grid;
	protected final Int2ObjectMap<TreeSet<ST>> supportOnCell;
	private final TreeMultimap<DistortedFace, ST> supportOnFace;
	protected TreeSet<ST> shapeFunctions;
	protected SparseMatrix systemMatrix;
	protected DenseVector rhs;
	Set<Integer> fixedNodes;
	private final double maxDiam;
	
	@Override
	public Set<Integer> getFixedNodeIndices()
	{
		return fixedNodes;
	}
	
	public DistortedGridSpace(final DistortedGrid grid)
	{
		supportOnCell = new Int2ObjectArrayMap<>();
		supportOnFace = TreeMultimap.create();
		this.grid = grid;
		fixedNodes = Sets.newConcurrentHashSet();
		maxDiam = getCells().stream()
		                    .mapToDouble(DistortedCell::diam)
		                    .max()
		                    .orElse(1);
	}
	
	public double getMaxDiam()
	{
		return maxDiam;
	}
	
	protected void addFunctionToCell(final ST function, final DistortedCell cell)
	{
		supportOnCell.putIfAbsent(cell.doneCode(), new TreeSet<>());
		supportOnCell.get(cell.doneCode())
		             .add(function);
	}
	
	protected void addFunctionToFace(final ST function, final DistortedFace face)
	{
		supportOnFace.put(face, function);
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
		return grid.getDimension();
	}
	
	@Override
	public List<DistortedCell> getCells()
	{
		return grid.getCells()
		           .asList();
	}
	
	@Override
	public Collection<ST> getShapeFunctions()
	{
		return shapeFunctions;
	}
	
	@Override
	public Map<Integer, ST> getShapeFunctionMap()
	{
		
		final Map<Integer, ST> functionNumbers = new TreeMap<>();
		for (final ST shapeFunction : shapeFunctions)
			functionNumbers.put(shapeFunction.getGlobalIndex(), shapeFunction);
		return functionNumbers;
	}
	
	@Override
	public List<DistortedFace> getFaces()
	{
		return grid.getFaces()
		           .asList();
	}
	
	@Override
	public Set<ST> getShapeFunctionsWithSupportOnCell(final DistortedCell cell)
	{
		return supportOnCell.get(cell.doneCode());
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
	
	@Override
	public Multimap<DistortedFace, ST> getFaceSupportMapping()
	{
		return supportOnFace;
	}
	
	public List<CoordinateVector> generateIsotropicPoints(final int resolution)
	{
		return grid.generateIsotropicPlotPoints(resolution);
	}
	
	public List<Tuple2<DistortedCell, CoordinateVector>> generateReferencePlotPoints(final int resolution)
	{
		return grid.generateReferencePlotPoints(resolution);
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
				final boolean isInSide = getCells().stream()
				                                   .anyMatch(c -> c.isInCellPrecise(pos));
				if (isInSide) return 1.;
				else return 0.;
			}
		};
	}
}
