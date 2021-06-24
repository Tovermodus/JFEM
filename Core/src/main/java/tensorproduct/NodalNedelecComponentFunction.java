package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.NodeFunctional;
import basic.ScalarFunction;
import basic.ScalarShapeFunction;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import linalg.Matrix;

import java.util.*;

public class NodalNedelecComponentFunction extends ScalarShapeFunction<TPCell, TPFace, TPEdge, NodalNedelecComponentFunction> implements Comparable<NodalNedelecComponentFunction> {
	
	private final Map<TPCell, List<RTBasisFunction1D>> cells;
	private final Set<TPFace> faces;
	private final LagrangeNodeFunctional nodeFunctional;
	private final int polynomialDegree;
	private final int lowDegreeDimension;
	private final int dimension;
	public NodalNedelecComponentFunction(TPCell supportCell, int polynomialDegree, int localIndex,
	                                     int lowDegreeDimension)
	{
		cells = new TreeMap<>();
		faces = new TreeSet<>();
		this.polynomialDegree = polynomialDegree;
		this.lowDegreeDimension = lowDegreeDimension;
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
				cell.cell1Ds.get(i), i != lowDegreeDimension));
		}
		return function1Ds;
	}
	private List<RTBasisFunction1D> generateBasisFunctionOnCell(TPCell cell,
	                                                                  CoordinateVector functionalPoint)
	{
		if(!cell.isInCell(functionalPoint))
			throw new IllegalArgumentException("functional point is not in cell");
		List<RTBasisFunction1D> function1Ds = new ArrayList<>();
		for (int i = 0; i < functionalPoint.getLength(); i++)
		{
			function1Ds.add(new RTBasisFunction1D(polynomialDegree, functionalPoint.at(i),
				cell.cell1Ds.get(i), i != lowDegreeDimension));
		}
		return function1Ds;
	}
	private int[] decomposeIndex(int dimension, int polynomialDegree, int localIndex)
	{
		int[] ret = new int[dimension];
		for (int i = 0; i < dimension; i++)
		{
			if(i != lowDegreeDimension)
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
	public List<RTBasisFunction1D> get1DFunctionsInCell(TPCell cell)
	{
		return cells.get(cell);
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
	public NodeFunctional<ScalarFunction, Double, CoordinateVector, Matrix> getNodeFunctional() {
		return nodeFunctional;
	}
	
	@Override
	public void addFace(TPFace face) {
		faces.add(face);
	}
	
	@Override
	public void addCell(TPCell cell) {
		if(!cells.containsKey(cell))
		{
			cells.put(cell, generateBasisFunctionOnCell(cell, nodeFunctional.getPoint()));
		}
	}
	
	
	@Override
	public boolean hasFastEvaluation() {
		return true;
	}
	
	@Override
	public Double value(CoordinateVector pos)
	{
		for(TPCell  c:cells.keySet())
		{
			if(c.isInCell(pos))
				return valueInCell(pos,c);
		}
		return 0.;
	}
	
	@Override
	public CoordinateVector gradient(CoordinateVector pos)
	{
		for(TPCell  c:cells.keySet())
		{
			if(c.isInCell(pos))
				return gradientInCell(pos,c);
		}
		return new CoordinateVector(pos.getLength());
	}
	@Override
	public double fastValue(CoordinateVector pos) {
		for(TPCell  c:cells.keySet())
		{
			if(c.isInCell(pos))
				return fastValueInCell(pos,c);
		}
		return 0.;
	}
	
	@Override
	public double[] fastGradient(CoordinateVector pos) {
		for(TPCell  c:cells.keySet())
		{
			if(c.isInCell(pos))
				return fastGradientInCell(pos,c);
		}
		return new double[pos.getLength()];
	}
	
	@Override
	public Double valueInCell(CoordinateVector pos, TPCell cell) {
		return fastValueInCell(pos,cell);
	}
	
	@Override
	public CoordinateVector gradientInCell(CoordinateVector pos, TPCell cell) {
		return CoordinateVector.fromValues(fastGradientInCell(pos,cell));
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
	public int compareTo(NodalNedelecComponentFunction o) {
		return CoordinateComparator.comp(nodeFunctional.getPoint(), o.nodeFunctional.getPoint());
	}
	
	@Override
	public String toString() {
		return "Cell: ".concat(", Node point: ").concat(nodeFunctional.getPoint().toString()).concat(", global Index: ").concat(getGlobalIndex()+"");
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(obj instanceof NodalNedelecComponentFunction)
			return CoordinateComparator.comp(nodeFunctional.getPoint(),
				((NodalNedelecComponentFunction) obj).nodeFunctional.getPoint()) == 0;
		return false;
	}
	
}
