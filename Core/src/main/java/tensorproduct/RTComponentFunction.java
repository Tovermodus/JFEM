package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.NodeFunctional;
import basic.PerformanceArguments;
import basic.ScalarShapeFunction;
import linalg.CoordinateComparator;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;

public class RTComponentFunction
	implements ScalarShapeFunction<TPCell, TPFace>,
	Comparable<RTComponentFunction>
{
	
	private final Map<TPCell, List<RTBasisFunction1D>> cells;
	private final Set<TPFace> faces;
	private final LagrangeNodeFunctional nodeFunctional;
	private final int polynomialDegree;
	private final int localIndex;
	private final int highDegreeDimension;
	private final int dimension;
	private int globalIndex;
	
	public RTComponentFunction(final TPCell supportCell, final int polynomialDegree, final int localIndex, final int highDegreeDimension)
	{
		cells = new TreeMap<>();
		faces = new TreeSet<>();
		this.polynomialDegree = polynomialDegree;
		this.localIndex = localIndex;
		this.highDegreeDimension = highDegreeDimension;
		dimension = supportCell.getDimension();
		final List<RTBasisFunction1D> supportCellFunctions = generateBasisFunctionOnCell(supportCell,
		                                                                                 localIndex);
		final CoordinateVector functionalPoint =
			CoordinateVector.fromValues(supportCellFunctions.stream()
			                                                .mapToDouble(RTBasisFunction1D::getDegreeOfFreedom)
			                                                .toArray());
		nodeFunctional = new LagrangeNodeFunctional(functionalPoint);
		cells.put(supportCell, supportCellFunctions);
		checkIfPointOnFace(functionalPoint, supportCell);
	}
	
	private void checkIfPointOnFace(final CoordinateVector functionalPoint, final TPCell cell)
	{
		
		for (final TPFace face : cell.getFaces())
		{
			if (faces.add(face))
			{
				if (face.isOnFace(functionalPoint))
				{
					for (final TPCell cellOfFace : face.getCells())
					{
						cells.put(cellOfFace,
						          generateBasisFunctionOnCell(cellOfFace, functionalPoint));
						checkIfPointOnFace(functionalPoint, cellOfFace);
					}
				}
			}
		}
	}
	
	private List<RTBasisFunction1D> generateBasisFunctionOnCell(final TPCell cell,
	                                                            final int localIndex)
	{
		final int[] decomposedLocalIndex = decomposeIndex(cell.getDimension(), polynomialDegree, localIndex);
		final List<RTBasisFunction1D> function1Ds = new ArrayList<>();
		for (int i = 0; i < getDomainDimension(); i++)
		{
			function1Ds.add(new RTBasisFunction1D(polynomialDegree, decomposedLocalIndex[i],
			                                      cell.getComponentCell(i), i == highDegreeDimension));
		}
		return function1Ds;
	}
	
	private List<RTBasisFunction1D> generateBasisFunctionOnCell(final TPCell cell,
	                                                            final CoordinateVector functionalPoint)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!cell.isInCell(functionalPoint))
				throw new IllegalArgumentException("functional point is not in cell");
		final List<RTBasisFunction1D> function1Ds = new ArrayList<>();
		for (int i = 0; i < dimension; i++)
		{
			function1Ds.add(new RTBasisFunction1D(polynomialDegree, functionalPoint.at(i),
			                                      cell.getComponentCell(i), i == highDegreeDimension));
		}
		return function1Ds;
	}
	
	private int[] decomposeIndex(final int dimension, final int polynomialDegree, int localIndex)
	{
		final int[] ret = new int[dimension];
		for (int i = 0; i < dimension; i++)
		{
			if (i == highDegreeDimension)
			{
				ret[i] = localIndex % (polynomialDegree + 2);
				localIndex = localIndex / (polynomialDegree + 2);
			} else
			{
				ret[i] = localIndex % (polynomialDegree + 1);
				localIndex = localIndex / (polynomialDegree + 1);
			}
		}
		return ret;
	}
	
	@Override
	public int getDomainDimension()
	{
		return dimension;
	}
	
	@Override
	public Set<TPCell> getCells()
	{
		
		return cells.keySet();
	}
	
	@Override
	public Set<TPFace> getFaces()
	{
		return faces;
	}
	
	@Override
	public NodeFunctional<Double, CoordinateVector, CoordinateMatrix> getNodeFunctional()
	{
		return nodeFunctional;
	}
	
	public void setGlobalIndex(final int index)
	{
		globalIndex = index;
	}
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	public double fastValueInCell(final CoordinateVector pos, final TPCell cell)
	{
		double ret = 1;
		if (cell == null)
			return ret;
		final List<? extends Function1D> function1Ds;
		if (cells.containsKey(cell))
		{
			function1Ds = cells.get(cell);
			for (int i = 0; i < pos.getLength(); i++)
			{
				ret *= function1Ds.get(i)
				                  .value(pos.at(i));
			}
			return ret;
		}
		return 0.;
	}
	
	public double[] fastGradientInCell(final CoordinateVector pos, final TPCell cell)
	{
		final double[] ret = new double[pos.getLength()];
		if (cell == null)
			return ret;
		final List<? extends Function1D> function1Ds;
		if (cells.containsKey(cell))
		{
			function1Ds = cells.get(cell);
			for (int i = 0; i < pos.getLength(); i++)
			{
				double component = 1;
				for (int j = 0; j < pos.getLength(); j++)
				{
					if (i == j)
						component *= function1Ds.get(j)
						                        .derivative(pos.at(j));
					else
						component *= function1Ds.get(j)
						                        .value(pos.at(j));
				}
				ret[i] = component;
			}
		}
		return ret;
	}
	
	@Override
	public Double valueInCell(final CoordinateVector pos, final TPCell cell)
	{
		return fastValueInCell(pos, cell);
	}
	
	@Override
	public CoordinateVector gradientInCell(final CoordinateVector pos, final TPCell cell)
	{
		return new CoordinateVector(fastGradientInCell(pos, cell));
	}
	
	@Override
	public int compareTo(final RTComponentFunction o)
	{
		return CoordinateComparator.comp(nodeFunctional.getPoint(), o.nodeFunctional.getPoint());
	}
	
	@Override
	public String toString()
	{
		return "Cell: ".concat(", Node point: ")
		               .concat(nodeFunctional.getPoint()
		                                     .toString())
		               .concat(", global Index: ")
		               .concat(getGlobalIndex() + "");
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (obj instanceof RTComponentFunction)
			return CoordinateComparator.comp(nodeFunctional.getPoint(),
			                                 ((RTComponentFunction) obj).nodeFunctional.getPoint()) == 0;
		return false;
	}
}
