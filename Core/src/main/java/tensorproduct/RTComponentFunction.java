package tensorproduct;

import basic.*;
import linalg.*;
import linalg.Vector;

import java.util.*;

public class RTComponentFunction implements FastEvaluatedScalarShapeFunction<TPCell, TPFace,TPEdge,RTComponentFunction>,
	Comparable<RTComponentFunction> {
	
	private Map<TPCell, List<RTBasisFunction1D>> cells;
	private Set<TPFace> faces;
	private final LagrangeNodeFunctional nodeFunctional;
	private final int polynomialDegree;
	private int localIndex;
	private final int highDegreeDimension;
	private final int dimension;
	private int globalIndex;
	
	public RTComponentFunction(TPCell supportCell, int polynomialDegree,int localIndex, int highDegreeDimension)
	{
		cells = new TreeMap<>();
		faces = new TreeSet<>();
		this.polynomialDegree = polynomialDegree;
		this.localIndex = localIndex;
		this.highDegreeDimension = highDegreeDimension;
		dimension = supportCell.getDimension();
		List<RTBasisFunction1D> supportCellFunctions = generateBasisFunctionOnCell(supportCell,
			localIndex);
		CoordinateVector functionalPoint =
			CoordinateVector.fromValues(supportCellFunctions.stream().mapToDouble(RTBasisFunction1D::getDegreeOfFreedom).toArray());
		nodeFunctional = new LagrangeNodeFunctional(functionalPoint);
		cells.put(supportCell, supportCellFunctions);
		checkIfPointOnFace(functionalPoint,supportCell);
	}
	private void checkIfPointOnFace(CoordinateVector functionalPoint, TPCell cell)
	{
		
		for(TPFace face: cell.faces)
		{
			if(faces.add(face))
			{
				if (face.isOnFace(functionalPoint))
				{
					for (TPCell cellOfFace : face.getCells())
					{
						cells.put(cellOfFace, generateBasisFunctionOnCell(cellOfFace, functionalPoint));
						checkIfPointOnFace(functionalPoint,cellOfFace);
					}
				}
			}
		}
	}
	private List<RTBasisFunction1D> generateBasisFunctionOnCell(TPCell cell,
	                                                                  int localIndex)
	{
		int[] decomposedLocalIndex = decomposeIndex(cell.getDimension(), polynomialDegree, localIndex);
		List<RTBasisFunction1D> function1Ds = new ArrayList<>();
		for (int i = 0; i < getDomainDimension(); i++)
		{
			function1Ds.add(new RTBasisFunction1D(polynomialDegree, decomposedLocalIndex[i],
				cell.cell1Ds.get(i), i == highDegreeDimension));
		}
		return function1Ds;
	}
	private List<RTBasisFunction1D> generateBasisFunctionOnCell(TPCell cell,
	                                                                  CoordinateVector functionalPoint)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		if(!cell.isInCell(functionalPoint))
			throw new IllegalArgumentException("functional point is not in cell");
		List<RTBasisFunction1D> function1Ds = new ArrayList<>();
		for (int i = 0; i < functionalPoint.getLength(); i++)
		{
			function1Ds.add(new RTBasisFunction1D(polynomialDegree, functionalPoint.at(i),
				cell.cell1Ds.get(i), i == highDegreeDimension));
		}
		return function1Ds;
	}
	private int[] decomposeIndex(int dimension, int polynomialDegree, int localIndex)
	{
		int[] ret = new int[dimension];
		for (int i = 0; i < dimension; i++)
		{
			if(i == highDegreeDimension)
			{
				ret[i] = localIndex % (polynomialDegree + 2);
				localIndex = localIndex / (polynomialDegree + 2);
			}
			else
			{
				ret[i] = localIndex % (polynomialDegree+1);
				localIndex = localIndex / (polynomialDegree+1);
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
	public Set<TPCell> getCells() {
		
		return cells.keySet();
	}
	
	@Override
	public Set<TPFace> getFaces() {
		return faces;
	}
	
	@Override
	public NodeFunctional<ScalarFunction, Double, CoordinateVector, CoordinateMatrix> getNodeFunctional() {
		return nodeFunctional;
	}
	
	public void setGlobalIndex(int index)
	{
		globalIndex = index;
	}
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	public void addFace(TPFace face) {
		
		if(true)
			throw new IllegalArgumentException();faces.add(face);
	}
	
	public void addCell(TPCell cell) {
		if(true)
			throw new IllegalArgumentException();
		if(!cells.containsKey(cell))
		{
			cells.put(cell, generateBasisFunctionOnCell(cell, nodeFunctional.getPoint()));
		}
	}
	
	@Override
	public double fastValueInCell(CoordinateVector pos, TPCell cell)
	{
		double ret = 1;
		if(cell == null)
			return ret;
		List<? extends Function1D> function1Ds;
		if(cells.containsKey(cell))
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
	public double[] fastGradientInCell(CoordinateVector pos, TPCell cell)
	{
		double[] ret = new double[pos.getLength()];
		if(cell == null)
			return ret;
		List<? extends Function1D> function1Ds;
		if(cells.containsKey(cell))
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
	public int compareTo(RTComponentFunction o) {
		return CoordinateComparator.comp(nodeFunctional.getPoint(), o.nodeFunctional.getPoint());
	}
	
	@Override
	public String toString() {
		return "Cell: ".concat(", Node point: ").concat(nodeFunctional.getPoint().toString()).concat(", global Index: ").concat(getGlobalIndex()+"");
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(obj instanceof RTComponentFunction)
			return CoordinateComparator.comp(nodeFunctional.getPoint(),
				((RTComponentFunction) obj).nodeFunctional.getPoint()) == 0;
		return false;
	}
	
}
