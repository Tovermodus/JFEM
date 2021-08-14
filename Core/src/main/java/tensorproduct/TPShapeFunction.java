package tensorproduct;

import basic.FastEvaluatedScalarShapeFunction;
import basic.LagrangeNodeFunctional;
import basic.ScalarShapeFunctionWithReferenceShapeFunction;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class TPShapeFunction implements ScalarShapeFunctionWithReferenceShapeFunction<TPCell, TPFace>,
	FastEvaluatedScalarShapeFunction<TPCell, TPFace>, Comparable<TPShapeFunction>
{
	private final int polynomialDegree;
	List<LagrangeBasisFunction1D> function1Ds;
	public final TPCell supportCell;
	LagrangeNodeFunctional nodeFunctional;
	private int globalIndex;
	
	public TPShapeFunction(final TPCell supportCell, final int polynomialDegree, final int localIndex)
	{
		this(supportCell, polynomialDegree, decomposeIndex(supportCell.getDimension(), polynomialDegree,
		                                                   localIndex));
	}
	
	public TPShapeFunction(final TPCell supportCell, final int polynomialDegree, final int[] localIndices)
	{
		this.polynomialDegree = polynomialDegree;
		this.supportCell = supportCell;
		function1Ds = generateBasisFunctionOnCell(supportCell, localIndices);
		final CoordinateVector functional_point =
			CoordinateVector.fromValues(function1Ds
				                            .stream()
				                            .mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom)
				                            .toArray());
		nodeFunctional = new LagrangeNodeFunctional(functional_point);
	}
	
	public TPShapeFunction(final TPCell supportCell, final int polynomialDegree, final CoordinateVector functionalPoint)
	{
		this.polynomialDegree = polynomialDegree;
		this.supportCell = supportCell;
		function1Ds = generateBasisFunctionOnCell(supportCell, functionalPoint);
		final CoordinateVector functional_point =
			CoordinateVector.fromValues(function1Ds
				                            .stream()
				                            .mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom)
				                            .toArray());
		nodeFunctional = new LagrangeNodeFunctional(functional_point);
	}
	
	public TPShapeFunction(final TPCell supportCell, final List<LagrangeBasisFunction1D> function1Ds)
	{
		this.polynomialDegree = function1Ds.get(0).getPolynomialDegree();
		this.function1Ds = function1Ds;
		this.supportCell = supportCell;
		final CoordinateVector functional_point =
			CoordinateVector.fromValues(function1Ds
				                            .stream()
				                            .mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom)
				                            .toArray());
		nodeFunctional = new LagrangeNodeFunctional(functional_point);
	}
	
	private List<LagrangeBasisFunction1D> generateBasisFunctionOnCell(final TPCell cell,
	                                                                  final int[] localIndices)
	{
		final List<LagrangeBasisFunction1D> function1Ds = new ArrayList<>();
		for (int i = 0; i < localIndices.length; i++)
		{
			function1Ds.add(new LagrangeBasisFunction1D(polynomialDegree, localIndices[i],
			                                            cell.getComponentCell(i)));
		}
		return function1Ds;
	}
	
	private List<LagrangeBasisFunction1D> generateBasisFunctionOnCell(final TPCell cell,
	                                                                  final CoordinateVector functionalPoint)
	{
		if (!cell.isInCell(functionalPoint))
			throw new IllegalArgumentException(
				"functional point is not in cell" + functionalPoint + " " + cell);
		final List<LagrangeBasisFunction1D> function1Ds = new ArrayList<>();
		for (int i = 0; i < functionalPoint.getLength(); i++)
		{
			function1Ds.add(new LagrangeBasisFunction1D(polynomialDegree, functionalPoint.at(i),
			                                            cell.getComponentCell(i)));
		}
		return function1Ds;
	}
	
	public static int functionsPerCell(final int polynomialDegree, final int dimension)
	{
		return (int) Math.pow(polynomialDegree + 1, dimension);
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
	
	@Override
	public Set<TPCell> getCells()
	{
		final Set<TPCell> ret = new HashSet<>();
		ret.add(supportCell);
		return ret;
	}
	
	@Override
	public Set<TPFace> getFaces()
	{
		return supportCell.getFaces();
	}
	
	@Override
	public LagrangeNodeFunctional getNodeFunctional()
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
	
	@Override
	public double fastValueInCell(final CoordinateVector pos, final TPCell cell)
	{
		if (cell == null)
			return 0.;
		if (cell.equals(supportCell))
		{
			return fastValue(pos);
		} else
			return 0.;
	}
	
	@Override
	public double[] fastGradientInCell(final CoordinateVector pos, final TPCell cell)
	{
		if (cell == null)
			return new double[pos.getLength()];
		if (cell.equals(supportCell))
			return fastGradient(pos);
		else
			return new double[pos.getLength()];
	}
	
	@Override
	public double fastValue(final CoordinateVector pos)
	{
		if (!supportCell.isInCell(pos))
			return 0;
		double ret = 1;
		for (int i = 0; i < pos.getLength(); i++)
		{
			ret *= function1Ds.get(i).value(pos.at(i));
		}
		return ret;
	}
	
	@Override
	public double[] fastGradient(final CoordinateVector pos)
	{
		if (!supportCell.isInCell(pos))
			return new double[pos.getLength()];
		final double[] ret = new double[pos.getLength()];
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
		return ret;
	}
	
	@Override
	public int hashCode()
	{
		return nodeFunctional.getPoint().hashCode();
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (obj instanceof TPShapeFunction)
			return 0 == compareTo((TPShapeFunction) obj);
		else
			return false;
	}
	
	@Override
	public String toString()
	{
		return "Cell: ".concat(supportCell.toString()).concat(", Node point: ").concat(
			nodeFunctional.getPoint().toString()).concat(", global Index: ").concat(getGlobalIndex() + "");
	}
	
	@Override
	public TPShapeFunction createReferenceShapeFunctionRelativeTo(final TPCell cell)
	{
		return new TPShapeFunction(cell.getReferenceCell(), polynomialDegree,
		                           function1Ds
			                           .stream()
			                           .mapToInt(LagrangeBasisFunction1D::getLocalFunctionNumber)
			                           .toArray());
	}
	
	@Override
	public TPShapeFunction createReferenceShapeFunctionRelativeTo(final TPFace face)
	{
		if (face.isNormalDownstream(supportCell.center()))
			return new TPShapeFunction(face.getReferenceFace().getNormalDownstreamCell(), polynomialDegree,
			                           function1Ds
				                           .stream()
				                           .mapToInt(LagrangeBasisFunction1D::getLocalFunctionNumber)
				                           .toArray());
		else
			return new TPShapeFunction(face.getReferenceFace().getNormalUpstreamCell(), polynomialDegree,
			                           function1Ds
				                           .stream()
				                           .mapToInt(LagrangeBasisFunction1D::getLocalFunctionNumber)
				                           .toArray());
	}
	
	@Override
	public int compareTo(final TPShapeFunction o)
	{
		if (polynomialDegree > o.polynomialDegree)
			return 1;
		if (polynomialDegree < o.polynomialDegree)
			return -1;
		if (CoordinateComparator.comp(nodeFunctional.getPoint().getEntries(),
		                              o.nodeFunctional.getPoint().getEntries()) == 0)
			return CoordinateComparator.comp(supportCell.center().getEntries(),
			                                 o.supportCell.center().getEntries());
		return CoordinateComparator.comp(nodeFunctional.getPoint().getEntries(),
		                                 o.nodeFunctional.getPoint().getEntries());
	}
}
