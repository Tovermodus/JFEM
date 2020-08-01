package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.NodeFunctional;
import basic.ScalarShapeFunction;
import basic.ShapeFunction;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import linalg.Matrix;
import linalg.Vector;
import org.jetbrains.annotations.NotNull;

import java.util.*;

public class TPShapeFunction extends ScalarShapeFunction<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,
	TPShapeFunction> implements Comparable<TPShapeFunction>
{
	List<LagrangeBasisFunction1D> function1Ds;
	TPCell<TPShapeFunction> supportCell;
	LagrangeNodeFunctional nodeFunctional;
	
	public TPShapeFunction(TPCell<TPShapeFunction> supportCell, int polynomialDegree, int localIndex)
	{
		this(supportCell, polynomialDegree, decomposeIndex(supportCell.getDimension(),polynomialDegree,
			localIndex));
	}
	public TPShapeFunction(TPCell<TPShapeFunction> supportCell, int polynomialDegree, int[] localIndices)
	{
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
		this.supportCell.addShapeFunction(this);
		for(TPFace<TPShapeFunction> f: getFaces())
			f.addShapeFunction(this);
	}
	public TPShapeFunction(TPCell<TPShapeFunction> supportCell, List<LagrangeBasisFunction1D> function1Ds)
	{
		this.function1Ds = function1Ds;
		this.supportCell = supportCell;
		CoordinateVector functional_point =
			CoordinateVector.fromValues(function1Ds.stream().mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom).toArray());
		nodeFunctional = new LagrangeNodeFunctional(functional_point);
		this.supportCell.addShapeFunction(this);
		for(TPFace<TPShapeFunction> f: getFaces())
			f.addShapeFunction(this);
	}
	private static int[] decomposeIndex(int dimension, int polynomialDegree, int localIndex)
	{
		int[] ret = new int[dimension];
		for (int i = 0; i < dimension; i++)
		{
			ret[i] = localIndex % polynomialDegree;
			localIndex = localIndex/polynomialDegree;
		}
		return ret;
	}
	@Override
	public Set<TPCell<TPShapeFunction>> getCells()
	{
		Set<TPCell<TPShapeFunction>> ret = new HashSet<TPCell<TPShapeFunction>>();
		ret.add(supportCell);
		return ret;
	}
	
	@Override
	public Set<TPFace<TPShapeFunction>> getFaces()
	{
		return supportCell.getFaces();
	}
	
	@Override
	public NodeFunctional getNodeFunctional()
	{
		return nodeFunctional;
	}
	
	@Override
	public void addFace(TPFace<TPShapeFunction> face)
	{
		supportCell.addFace(face);
		face.addShapeFunction(this);
	}
	
	@Override
	public void addCell(TPCell<TPShapeFunction> cell)
	{
		if(supportCell != cell)
		{
			if(supportCell != null)
				throw new UnsupportedOperationException();
			supportCell = cell;
			cell.addShapeFunction(this);
		}
	}
	
	@Override
	public Double valueInCell(CoordinateVector pos, TPCell<TPShapeFunction> cell)
	{
		if(cell == null)
			return 0.;
		if(cell.equals(supportCell))
			return value(pos);
		else
			return 0.;
	}
	
	@Override
	public Vector gradientInCell(CoordinateVector pos, TPCell<TPShapeFunction> cell)
	{
		if(cell == null)
			return gradient(pos).mul(0);
		if(cell.equals(supportCell))
			return gradient(pos);
		else
			return gradient(pos).mul(0);
	}
	
	@Override
	public Matrix hessianInCell(CoordinateVector pos, TPCell<TPShapeFunction> cell)
	{
		if(cell == null)
			return hessian(pos).mul(0);
		if(cell.equals(supportCell))
			return hessian(pos);
		else
			return hessian(pos).mul(0);
	}
	
	@Override
	public Map<Integer, Double> prolongate(Set<TPShapeFunction> refinedFunctions)
	{
		Map<Integer, Double> ret = new HashMap<>();
		for(TPShapeFunction shapeFunction:refinedFunctions)
		{
			ret.put(shapeFunction.getGlobalIndex(), shapeFunction.getNodeFunctional().evaluate(this));
		}
		return ret;
	}
	
	@Override
	public Double value(CoordinateVector pos)
	{
		double ret = 1;
		for(int i = 0; i < pos.getLength(); i++)
		{
			ret *= function1Ds.get(i).value(pos.at(i));
		}
		return ret;
	}
	
	@Override
	public Vector gradient(CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(pos.getLength());
		for (int i = 0; i < pos.getLength(); i++)
		{
			double comp = 1;
			for(int j = 0; j < pos.getLength(); j++)
			{
				if(i == j)
					comp *= function1Ds.get(j).derivative(pos.at(j));
				else
					comp *= function1Ds.get(j).value(pos.at(j));
			}
			ret.set(comp,i);
		}
		return ret;
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
			double comp = 1;
			for(int j = 0; j < pos.getLength(); j++)
			{
				if(i == j)
					comp *= function1Ds.get(i).derivative(pos.at(i));
				else
					comp *= function1Ds.get(i).value(pos.at(i));
				System.out.println("comp"+comp+" i "+i);
			}
			ret[i] = comp;
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
	//	LagrangeBasisFunction1D xFunction;
//	LagrangeBasisFunction1D yFunction;
//	private CoordinateVector functional_point;
//	TPCell supportCell;
//	Set<TPFace> faces;
//	LagrangeNodeFunctional nodeFunctional;
//	int globalIndex;
//	public TPShapeFunction(LagrangeBasisFunction1D xFunction, LagrangeBasisFunction1D yFunction,TPCell cell)
//	{
//		this.xFunction = xFunction;
//		this.yFunction = yFunction;
//		functional_point = CoordinateVector.vectorFromValues(xFunction.degreeOfFreedom, yFunction.degreeOfFreedom);
//		cell.getShapeFunctions().add(this);
//		for(TPFace f: cell.getFaces())
//			f.getShapeFunctions().add(this);
//		nodeFunctional = new LagrangeNodeFunctional(functional_point);
//	}
//
//	@Override
//	public Double value(CoordinateVector pos)
//	{
//		return xFunction.value(pos.x())*yFunction.value(pos.y());
//	}
//
//	@Override
//	public Vector derivative(CoordinateVector pos)
//	{
//		return .vectorFromValues(xFunction.derivative(pos.x())*yFunction.value(pos.y()),
//			xFunction.value(pos.x())*yFunction.derivative(pos.y()));
//	}
//
//	@Override
//	public void setGlobalIndex(int index)
//	{
//		globalIndex = index;
//	}
//
//	@Override
//	public List<TPCell> getCells()
//	{
//		Set<TPCell> ret = new HashSet<>();
//		ret.add(supportCell);
//		return ret;
//	}
//
//	@Override
//	public int getGlobalIndex()
//	{
//		return globalIndex;
//	}
//
//	@Override
//	public ScalarNodeFunctional getNodeFunctional()
//	{
//		return nodeFunctional;
//	}
//
//	@Override
//	public double valueInCell(DoubleTensor pos, TPCell cell)
//	{
//		if(supportCell == cell)
//			return value(pos);
//		else
//			return 0;
//	}
//	@Override
//	public DoubleTensor derivativeInCell(DoubleTensor pos, TPCell cell)
//	{
//		if(supportCell == cell)
//			return derivative(pos);
//		else
//			return new DoubleTensor(pos.size());
//	}
//
//	@Override
//	public Map<Integer, Double> prolongate(List<ScalarShapeFunction<TPCell, TPFace>> refinedFunctions)
//	{
//	}
//
//
//	@Override
//	public int compareTo(@NotNull TPShapeFunction o)
//	{
//
//	}
}
