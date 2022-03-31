package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.ScalarShapeFunction;
import basic.ScalarShapeFunctionWithReferenceShapeFunction;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

public class TPShapeFunction
	implements ScalarShapeFunctionWithReferenceShapeFunction<TPCell, TPFace>,
	ScalarShapeFunction<TPCell, TPFace>, Comparable<TPShapeFunction>
{
	private final int polynomialDegree;
	LagrangeBasisFunction1D[] function1Ds;
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
			CoordinateVector.fromValues(Arrays.stream(function1Ds)
			                                  .mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom)
			                                  .toArray());
		nodeFunctional = new LagrangeNodeFunctional(functional_point);
	}
	
	public TPShapeFunction(final TPCell supportCell,
	                       final int polynomialDegree,
	                       final CoordinateVector functionalPoint)
	{
		this.polynomialDegree = polynomialDegree;
		this.supportCell = supportCell;
		function1Ds = generateBasisFunctionOnCell(supportCell, functionalPoint);
		final CoordinateVector functional_point =
			CoordinateVector.fromValues(Arrays.stream(function1Ds)
			                                  .mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom)
			                                  .toArray());
		nodeFunctional = new LagrangeNodeFunctional(functional_point);
	}
	
	public TPShapeFunction(final TPCell supportCell, final List<LagrangeBasisFunction1D> function1Ds)
	{
		this.polynomialDegree = function1Ds.get(0)
		                                   .getPolynomialDegree();
		this.function1Ds = function1Ds.toArray(LagrangeBasisFunction1D[]::new);
		this.supportCell = supportCell;
		final CoordinateVector functional_point =
			CoordinateVector.fromValues(function1Ds
				                            .stream()
				                            .mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom)
				                            .toArray());
		nodeFunctional = new LagrangeNodeFunctional(functional_point);
	}
	
	static Int2ObjectMap<Int2ObjectMap<TPShapeFunction>> polynomialDegreeLocalIndexStore;
	
	public static TPShapeFunction ReferenceFunction(final int dimension,
	                                                final int polynomialDegree,
	                                                final int localIndex)
	{
		if (polynomialDegreeLocalIndexStore == null)
			polynomialDegreeLocalIndexStore = new Int2ObjectArrayMap<>();
		if (polynomialDegreeLocalIndexStore.get(polynomialDegree) == null)
			polynomialDegreeLocalIndexStore.put(polynomialDegree, new Int2ObjectArrayMap<>());
		final Int2ObjectMap<TPShapeFunction> thisStore = polynomialDegreeLocalIndexStore.get(polynomialDegree);
		if (!thisStore.containsKey(localIndex))
			thisStore.put(localIndex, new TPShapeFunction(TPCell.unitHyperCube(dimension),
			                                              polynomialDegree, localIndex));
		return thisStore.get(localIndex);
	}
	
	public static int getIndexFromReferenceNodePoint(final int polynomialDegree,
	                                                 final CoordinateVector referenceNodePoint)
	{
		final int dimension = referenceNodePoint.getLength();
		final int[] localIndices = new int[dimension];
		for (int i = 0; i < dimension; i++)
		{
			localIndices[i] =
				LagrangeBasisFunction1D.getLocalFunctionNumberOnReferenceCell(polynomialDegree,
				                                                              referenceNodePoint.at(i));
		}
		final int localIndex = composeIndex(dimension, polynomialDegree, localIndices);
		return localIndex;
	}
	
	private LagrangeBasisFunction1D[] generateBasisFunctionOnCell(final TPCell cell,
	                                                              final int[] localIndices)
	{
		final List<LagrangeBasisFunction1D> function1Ds = new ArrayList<>();
		for (int i = 0; i < localIndices.length; i++)
		{
			function1Ds.add(new LagrangeBasisFunction1D(polynomialDegree, localIndices[i],
			                                            cell.getComponentCell(i)));
		}
		return function1Ds.toArray(LagrangeBasisFunction1D[]::new);
	}
	
	private LagrangeBasisFunction1D[] generateBasisFunctionOnCell(final TPCell cell,
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
		return function1Ds.toArray(LagrangeBasisFunction1D[]::new);
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
	
	private static int composeIndex(final int dimension, final int polynomialDegree, final int[] localIndices)
	{
		int ret = 0;
		for (int i = dimension - 1; i >= 0; i--)
		{
			ret *= polynomialDegree + 1;
			ret += localIndices[i];
		}
		return ret;
	}
	
	@Override
	public List<TPCell> getCells()
	{
		return List.of(supportCell);
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
	
	public double fastValue(final CoordinateVector pos)
	{
		double ret = 1;
		for (int i = 0; i < pos.getLength(); i++)
		{
			ret *= function1Ds[i]
				.value(pos.at(i));
		}
		return ret;
	}
	
	public double fastValueOnReferenceCell(final CoordinateVector pos)
	{
		double ret = 1;
		for (int i = 0; i < pos.getLength(); i++)
		{
			ret *= function1Ds[i]
				.valueOnReferenceCell(pos.at(i));
		}
		return ret;
	}
	
	public double[] fastGradient(final CoordinateVector pos)
	{
		final int d = pos.getLength();
		final double[] ret = new double[d];
		for (int i = 0; i < d; i++)
		{
			double component = 1;
			for (int j = 0; j < d; j++)
			{
				if (i == j)
					component *= function1Ds[j]
						.derivative(pos.at(j));
				else
					component *= function1Ds[j]
						.value(pos.at(j));
			}
			ret[i] = component;
		}
		return ret;
	}
	
	public double[] fastGradientOnReferenceCell(final CoordinateVector pos)
	{
		final int d = pos.getLength();
		final double[] ret = new double[d];
		for (int i = 0; i < d; i++)
		{
			double component = 1;
			for (int j = 0; j < d; j++)
			{
				if (i == j)
					component *= function1Ds[j]
						.derivativeOnReferenceCell(pos.at(j));
				else
					component *= function1Ds[j]
						.valueOnReferenceCell(pos.at(j));
			}
			ret[i] = component;
		}
		return ret;
	}
	
	@Override
	public Double value(final CoordinateVector pos)
	{
		return fastValue(pos);
	}
	
	@Override
	public CoordinateVector gradient(final CoordinateVector pos)
	{
		return new CoordinateVector(fastGradient(pos));
	}
	
	@Override
	public Double valueInCell(final CoordinateVector pos, final TPCell cell)
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
	public CoordinateVector gradientInCell(final CoordinateVector pos, final TPCell cell)
	{
		if (cell == null)
			return defaultGradient();
		if (cell.equals(supportCell))
		{
			return new CoordinateVector(fastGradient(pos));
		} else
			return defaultGradient();
	}
	
	@Override
	public int hashCode()
	{
		return nodeFunctional.getPoint()
		                     .hashCode();
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
		return "Cell: ".concat(supportCell.toString())
		               .concat(", Node point: ")
		               .concat(
			               nodeFunctional.getPoint()
			                             .toString())
		               .concat(", global Index: ")
		               .concat(getGlobalIndex() + "");
	}
	
	@Override
	public TPShapeFunction createReferenceShapeFunctionRelativeTo(final TPCell cell)
	{
		return new TPShapeFunction(cell.getReferenceCell(), polynomialDegree, Arrays.stream(function1Ds)
		                                                                            .mapToInt(
			                                                                            LagrangeBasisFunction1D::getLocalFunctionNumber)
		                                                                            .toArray());
	}
	
	@Override
	public TPShapeFunction createReferenceShapeFunctionRelativeTo(final TPFace face)
	{
		if (face.isNormalDownstream(supportCell.center()))
			return new TPShapeFunction(face.getReferenceFace()
			                               .getNormalDownstreamCell(),
			                           polynomialDegree,
			                           Arrays.stream(function1Ds)
			                                 .mapToInt(LagrangeBasisFunction1D::getLocalFunctionNumber)
			                                 .toArray());
		else
			return new TPShapeFunction(face.getReferenceFace()
			                               .getNormalUpstreamCell(),
			                           polynomialDegree,
			                           Arrays.stream(function1Ds)
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
		if (CoordinateComparator.comp(nodeFunctional.getPoint()
		                                            .getEntries(),
		                              o.nodeFunctional.getPoint()
		                                              .getEntries()) == 0)
			return CoordinateComparator.comp(supportCell.center()
			                                            .getEntries(),
			                                 o.supportCell.center()
			                                              .getEntries());
		return CoordinateComparator.comp(nodeFunctional.getPoint()
		                                               .getEntries(),
		                                 o.nodeFunctional.getPoint()
		                                                 .getEntries());
	}
}
