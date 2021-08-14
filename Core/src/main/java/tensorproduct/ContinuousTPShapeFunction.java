package tensorproduct;

import basic.FastEvaluatedScalarShapeFunction;
import basic.LagrangeNodeFunctional;
import basic.NodeFunctional;
import basic.ScalarShapeFunctionWithReferenceShapeFunction;
import linalg.CoordinateComparator;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;
import java.util.stream.IntStream;

public class ContinuousTPShapeFunction implements ScalarShapeFunctionWithReferenceShapeFunction<TPCell, TPFace>, FastEvaluatedScalarShapeFunction<TPCell, TPFace>,
	Comparable<ContinuousTPShapeFunction>
{
	
	private final Map<TPCell, List<LagrangeBasisFunction1D>> cells;
	private final Set<TPFace> faces;
	private final LagrangeNodeFunctional nodeFunctional;
	private final int polynomialDegree;
	private int globalIndex;
	
	public ContinuousTPShapeFunction(final TPCell supportCell, final int polynomialDegree, final int localIndex)
	{
		cells = new TreeMap<>();
		faces = new TreeSet<>();
		this.polynomialDegree = polynomialDegree;
		final List<LagrangeBasisFunction1D> supportCellFunctions = generateBasisFunctionOnCell(supportCell,
                                                                                                       localIndex);
		final CoordinateVector functionalPoint =
			CoordinateVector.fromValues(supportCellFunctions
				                            .stream()
				                            .mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom)
				                            .toArray());
		nodeFunctional = new LagrangeNodeFunctional(functionalPoint);
		cells.put(supportCell, supportCellFunctions);
		checkIfPointOnFace(functionalPoint, supportCell);
	}
	
	public ContinuousTPShapeFunction(final TPCell supportCell, final int polynomialDegree, final CoordinateVector functionalPoint)
	{
		cells = new TreeMap<>();
		faces = new TreeSet<>();
		this.polynomialDegree = polynomialDegree;
		final List<LagrangeBasisFunction1D> supportCellFunctions = generateBasisFunctionOnCell(supportCell,
                                                                                                       functionalPoint);
		nodeFunctional = new LagrangeNodeFunctional(functionalPoint);
		cells.put(supportCell, supportCellFunctions);
		checkIfPointOnFace(functionalPoint, supportCell);
	}
	
	private static int[] decomposeIndex(final int dimension, final int polynomialDegree, int localIndex)
	{
		final int[] ret = new int[dimension];
		for (int i = 0; i < dimension; i++)
		{
			ret[i] = localIndex % (polynomialDegree + 1);
			localIndex = localIndex / (polynomialDegree + 1);
		}
		return ret;
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
	
	private List<LagrangeBasisFunction1D> generateBasisFunctionOnCell(final TPCell cell,
                                                                          final int localIndex)
	{
		final int[] decomposedLocalIndex = decomposeIndex(cell.getDimension(), polynomialDegree, localIndex);
		final List<LagrangeBasisFunction1D> function1Ds = new ArrayList<>();
		for (int i = 0; i < decomposedLocalIndex.length; i++)
		{
			function1Ds.add(new LagrangeBasisFunction1D(polynomialDegree, decomposedLocalIndex[i],
			                                            cell.getComponentCell(i)));
		}
		return function1Ds;
	}
	
	private List<LagrangeBasisFunction1D> generateBasisFunctionOnCell(final TPCell cell,
                                                                          final CoordinateVector functionalPoint)
	{
		if (!cell.isInCell(functionalPoint))
			throw new IllegalArgumentException("functional point is not in cell");
		final List<LagrangeBasisFunction1D> function1Ds = new ArrayList<>();
		for (int i = 0; i < functionalPoint.getLength(); i++)
		{
			function1Ds.add(new LagrangeBasisFunction1D(polynomialDegree, functionalPoint.at(i),
			                                            cell.getComponentCell(i)));
		}
		return function1Ds;
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
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	public void setGlobalIndex(final int index)
	{
		globalIndex = index;
	}
	
	@Override
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
				ret *= function1Ds.get(i).value(pos.at(i));
			}
			return ret;
		}
		return 0.;
	}
	
	@Override
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
						component *= function1Ds.get(j).derivative(pos.at(j));
					else
						component *= function1Ds.get(j).value(pos.at(j));
				}
				ret[i] = component;
			}
		}
		return ret;
	}
	
	@Override
	public String toString()
	{
		return "Cells: ".concat(getCells().toString()).concat(", Node point: ").concat(
			nodeFunctional.getPoint().toString()).concat(", " +
				                                             "global " +
				                                             "Index: ").concat(getGlobalIndex() + "");
	}
	
	@Override
	public int hashCode()
	{
		return nodeFunctional.getPoint().hashCode();
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (obj instanceof ContinuousTPShapeFunction)
			return CoordinateComparator.comp(nodeFunctional.getPoint(),
			                                 ((ContinuousTPShapeFunction) obj).nodeFunctional.getPoint()) ==
				0;
		return false;
	}
	
	@Override
        public int compareTo(final ContinuousTPShapeFunction o)
	{
		if (polynomialDegree > o.polynomialDegree)
			return 1;
		if (polynomialDegree < o.polynomialDegree)
			return -1;
		return CoordinateComparator.comp(nodeFunctional.getPoint().getEntries(),
		                                 o.nodeFunctional.getPoint().getEntries());
	}
	
	@Override
	public ContinuousTPShapeFunction createReferenceShapeFunctionRelativeTo(final TPCell cell)
	{
		final List<LagrangeBasisFunction1D> functions = cells.get(cell);
		final CoordinateVector functionalPoint =
			CoordinateVector.fromValues(
				IntStream.range(0, getDomainDimension())
				         .mapToDouble(i -> functions
					         .get(i)
					         .getCell()
					         .positionOnReferenceCell(nodeFunctional.getPoint().at(i)))
				         .toArray());
		return new ContinuousTPShapeFunction(cell.getReferenceCell(), polynomialDegree, functionalPoint);
	}
	
	@Override
	public ContinuousTPShapeFunction createReferenceShapeFunctionRelativeTo(final TPFace face)
	{
		final boolean functionalPointIsDownStreamOfFace =
			face.getNormalUpstreamCell() == null || face.isNormalDownstream(nodeFunctional.getPoint());
		final TPCell supportCell = functionalPointIsDownStreamOfFace ? face.getNormalDownstreamCell() :
                                           face.getNormalUpstreamCell();
		final TPCell referenceCell =
			functionalPointIsDownStreamOfFace ? face.getReferenceFace().getNormalDownstreamCell() :
			face.getReferenceFace().getNormalUpstreamCell();
		final List<LagrangeBasisFunction1D> functions = cells.get(supportCell);
		final CoordinateVector functionalPoint =
			CoordinateVector.fromValues(
				IntStream.range(0, getDomainDimension())
				         .mapToDouble(i -> functions
					         .get(i)
					         .getCell()
					         .positionOnReferenceCell(nodeFunctional.getPoint().at(i)))
				         .toArray());
		if (!functionalPointIsDownStreamOfFace)
			functionalPoint.add(-1, face.flatDimension);
		return new ContinuousTPShapeFunction(referenceCell, polynomialDegree, functionalPoint);
	}
}
