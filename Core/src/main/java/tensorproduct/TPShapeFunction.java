package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.NodeFunctional;
import basic.ScalarFunction;
import basic.ScalarShapeFunctionWithReferenceShapeFunction;
import linalg.CoordinateComparator;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class TPShapeFunction implements ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,TPEdge,
	TPShapeFunction>, Comparable<TPShapeFunction>
{
	private final int polynomialDegree;
	List<LagrangeBasisFunction1D> function1Ds;
	TPCell supportCell;
	LagrangeNodeFunctional nodeFunctional;
	private int globalIndex;
	
	public TPShapeFunction(TPCell supportCell, int polynomialDegree, int localIndex)
	{
		this(supportCell, polynomialDegree, decomposeIndex(supportCell.getDimension(),polynomialDegree,
			localIndex));
	}
	public TPShapeFunction(TPCell supportCell, int polynomialDegree, int[] localIndices)
	{
		this.polynomialDegree = polynomialDegree;
		function1Ds = new ArrayList<>();
		this.supportCell = supportCell;
		for (int i = 0; i < localIndices.length; i++)
		{
			function1Ds.add(new LagrangeBasisFunction1D(polynomialDegree, localIndices[i],
				this.supportCell.cell1Ds.get(i)));
		}
		CoordinateVector functional_point =
			CoordinateVector.fromValues(function1Ds.stream().mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom).toArray());
		nodeFunctional = new LagrangeNodeFunctional(functional_point);
	}
	public TPShapeFunction(TPCell supportCell, List<LagrangeBasisFunction1D> function1Ds)
	{
		this.polynomialDegree = function1Ds.get(0).getPolynomialDegree();
		this.function1Ds = function1Ds;
		this.supportCell = supportCell;
		CoordinateVector functional_point =
			CoordinateVector.fromValues(function1Ds.stream().mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom).toArray());
		nodeFunctional = new LagrangeNodeFunctional(functional_point);
	}
	private static int[] decomposeIndex(int dimension, int polynomialDegree, int localIndex)
	{
		int[] ret = new int[dimension];
		for (int i = 0; i < dimension; i++)
		{
			ret[i] = localIndex % (polynomialDegree+1);
			localIndex = localIndex/(polynomialDegree+1);
		}
		return ret;
	}
	@Override
	public Set<TPCell> getCells()
	{
		Set<TPCell> ret = new HashSet<>();
		ret.add(supportCell);
		return ret;
	}
	
	@Override
	public Set<TPFace> getFaces()
	{
		return supportCell.getFaces();
	}
	
	@Override
	public NodeFunctional<ScalarFunction, Double, CoordinateVector, CoordinateMatrix> getNodeFunctional()
	{
		return nodeFunctional;
	}
	
	@Override
	public void setGlobalIndex(int index)
	{
		globalIndex = index;
	}
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	@Override
	public void addFace(TPFace face)
	{
		throw new UnsupportedOperationException("TPShapefunctions do not live on faces");
	}
	
	@Override
	public void addCell(TPCell cell)
	{
		if(supportCell != cell)
		{
			if(supportCell != null)
				throw new UnsupportedOperationException();
			supportCell = cell;
		}
	}
	
	@Override
	public Double valueInCell(CoordinateVector pos, TPCell cell)
	{
		if(cell == null)
			return 0.;
		if(cell.equals(supportCell))
			return value(pos);
		else
			return 0.;
	}
	
	@Override
	public CoordinateVector gradientInCell(CoordinateVector pos, TPCell cell)
	{
		if(cell == null)
			return gradient(pos).mul(0);
		if(cell.equals(supportCell))
			return gradient(pos);
		else
			return gradient(pos).mul(0);
	}
	
	@Override
	public CoordinateMatrix hessianInCell(CoordinateVector pos, TPCell cell)
	{
		if(cell == null)
			return hessian(pos).mul(0);
		if(cell.equals(supportCell))
			return hessian(pos);
		else
			return hessian(pos).mul(0);
	}
	@Override
	public double fastValueInCell(CoordinateVector pos, TPCell cell)
	{
		if(cell == null)
			return 0.;
		if(cell.equals(supportCell))
			return fastValue(pos);
		else
			return 0.;
	}
	
	@Override
	public double[] fastGradientInCell(CoordinateVector pos, TPCell cell)
	{
		if(cell == null)
			return new double[pos.getLength()];
		if(cell.equals(supportCell))
			return fastGradient(pos);
		else
			return new double[pos.getLength()];
	}
	
	@Override
	public double[][] fastHessianInCell(CoordinateVector pos, TPCell cell)
	{
		if(cell == null)
			return new double[pos.getLength()][pos.getLength()];
		if(cell.equals(supportCell))
			return fastHessian(pos);
		else
			return new double[pos.getLength()][pos.getLength()];
	}
	
	
	@Override
	public Double value(CoordinateVector pos)
	{
		return fastValue(pos);
	}
	
	@Override
	public CoordinateVector gradient(CoordinateVector pos)
	{
		return CoordinateVector.fromValues(fastGradient(pos));
	}
	
	@Override
	public boolean hasFastEvaluation()
	{
		return true;
	}
	
	@Override
	public double fastValue(CoordinateVector pos)
	{
		double ret = 1;
		for(int i = 0; i < pos.getLength(); i++)
		{
			ret *= function1Ds.get(i).value(pos.at(i));
		}
		return ret;
	}
	
	@Override
	public double[] fastGradient(CoordinateVector pos)
	{
		double[] ret = new double[pos.getLength()];
		for (int i = 0; i < pos.getLength(); i++)
		{
			double component = 1;
			for(int j = 0; j < pos.getLength(); j++)
			{
				if(i == j)
					component *= function1Ds.get(j).derivative(pos.at(j));
				else
					component *= function1Ds.get(j).value(pos.at(j));
			}
			ret[i] = component;
		}
		return ret;
	}
	
	@Override
	public int compareTo(@NotNull TPShapeFunction o)
	{
		if(CoordinateComparator.comp(nodeFunctional.getPoint().getEntries(),
			o.nodeFunctional.getPoint().getEntries())==0)
			return CoordinateComparator.comp(supportCell.center().getEntries(),
				o.supportCell.center().getEntries());
		return CoordinateComparator.comp(nodeFunctional.getPoint().getEntries(),
			o.nodeFunctional.getPoint().getEntries());
	}
	
	@Override
	public String toString()
	{
		return "Cell: ".concat(supportCell.toString()).concat(", Node point: ").concat(nodeFunctional.getPoint().toString()).concat(", global Index: ").concat(getGlobalIndex()+"");
	}
	
	@Override
	public TPShapeFunction getReferenceShapeFunctionRelativeTo(TPCell cell)
	{
		return new TPShapeFunction(cell, polynomialDegree,
			function1Ds.stream().mapToInt(LagrangeBasisFunction1D::getLocalFunctionNumber).toArray());
	}
	
}
