package mixed;

import basic.*;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Tensor;
import org.jetbrains.annotations.NotNull;
import tensorproduct.*;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

public class NedelecNodeFuncional implements NodeFunctional<VectorFunction, CoordinateVector, CoordinateMatrix,
	Tensor>, Comparable<NedelecNodeFuncional>
{
	public static final int EDGETYPE = 0;
	public static final int FACETYPE = 1;
	public static final int CELLTYPE = 2;
	private final TPEdge e;
	private final TPFace f;
	
	private final TPCell c;
	Collection<TPCell> cells;
	List<LagrangeBasisFunction1D> testFunctions;
	int testComponent;
	
	public NedelecNodeFuncional(int polynomialDegree, TPEdge e, int number)
	{
		this.e = e;
		this.f = null;
		this.c = null;
		this.cells = e.getCells();
		testFunctions = List.of(new LagrangeBasisFunction1D(polynomialDegree - 1, number, e.getCell()));
		testComponent = 0;
	}
	
	public NedelecNodeFuncional(int polynomialDegree, TPFace f, int number)
	{
		this.e = null;
		this.f = f;
		this.c = null;
		testComponent = number / ((polynomialDegree - 1) * polynomialDegree);
		int numberInComponent = number % ((polynomialDegree - 1) * polynomialDegree);
		int testFunctionIndex1 = numberInComponent / (polynomialDegree - 1); //from 0 to polynomialDegree
		int testFunctionIndex2 = numberInComponent % (polynomialDegree - 1);// from 0 to polynomialDegree- 1
		this.cells = f.getCells();
		if (testComponent == 0)
			testFunctions = List.of(
				new LagrangeBasisFunction1D(polynomialDegree - 2, testFunctionIndex2,
					f.getCell1Ds().get(0)),
				new LagrangeBasisFunction1D(polynomialDegree - 1, testFunctionIndex1,
					f.getCell1Ds().get(1)));
		else
			testFunctions = List.of(
				new LagrangeBasisFunction1D(polynomialDegree - 1, testFunctionIndex2,
					f.getCell1Ds().get(0)),
				new LagrangeBasisFunction1D(polynomialDegree - 2, testFunctionIndex1,
					f.getCell1Ds().get(1)));
		
	}
	
	public NedelecNodeFuncional(int polynomialDegree, TPCell c, int number)
	{
		this.e = null;
		this.f = null;
		this.c = c;
		testComponent = number / ((polynomialDegree - 1) * (polynomialDegree - 1) * polynomialDegree);
		int numberInComponent = number % ((polynomialDegree - 1) * (polynomialDegree - 1) * polynomialDegree);
		int testFunctionIndex1 = numberInComponent / ((polynomialDegree - 1) * (polynomialDegree - 1)); //from 0 to
		// polynomialDegree
		numberInComponent = numberInComponent % ((polynomialDegree - 1) * (polynomialDegree - 1));
		int testFunctionIndex2 = numberInComponent % (polynomialDegree - 1);// from 0 to polynomialDegree- 1
		int testFunctionIndex3 = numberInComponent / (polynomialDegree - 1);// from 0 to polynomialDegree- 1
		this.cells = List.of(c);
		if (testComponent == 0)
			testFunctions = List.of(
				new LagrangeBasisFunction1D(polynomialDegree - 1, testFunctionIndex1,
					c.getCell1Ds().get(0)),
				new LagrangeBasisFunction1D(polynomialDegree - 2, testFunctionIndex2,
					c.getCell1Ds().get(1)),
				new LagrangeBasisFunction1D(polynomialDegree - 2, testFunctionIndex3,
					c.getCell1Ds().get(2)));
		if (testComponent == 1)
			testFunctions = List.of(
				new LagrangeBasisFunction1D(polynomialDegree - 2, testFunctionIndex2,
					c.getCell1Ds().get(0)),
				new LagrangeBasisFunction1D(polynomialDegree - 1, testFunctionIndex1,
					c.getCell1Ds().get(1)),
				new LagrangeBasisFunction1D(polynomialDegree - 2, testFunctionIndex3,
					c.getCell1Ds().get(2)));
		if (testComponent == 2)
			testFunctions = List.of(
				new LagrangeBasisFunction1D(polynomialDegree - 2, testFunctionIndex3,
					c.getCell1Ds().get(0)),
				new LagrangeBasisFunction1D(polynomialDegree - 2, testFunctionIndex2,
					c.getCell1Ds().get(1)),
				new LagrangeBasisFunction1D(polynomialDegree - 1, testFunctionIndex3,
					c.getCell1Ds().get(2)));
	}
	public Collection<TPCell> getCells()
	{
		return cells;
	}
	
	public TPEdge getE()
	{
		
		if(getType() != EDGETYPE)
			throw new IllegalStateException("not an edge functional");
		return e;
	}
	
	public TPFace getF()
	{
		if(getType() != FACETYPE)
			throw new IllegalStateException("not a face functional");
		return f;
	}
	
	public TPCell getC()
	{
		if(getType() != CELLTYPE)
			throw new IllegalStateException("not a cell functional");
		return c;
	}
	
	public int getType()
	{
		if (e != null)
			return EDGETYPE;
		else if (f != null)
			return FACETYPE;
		else
			return CELLTYPE;
	}
	
	@Override
	public double evaluate(VectorFunction func)
	{
		if (e != null)
		{
			return TPEdgeIntegral.integrateNonTensorProduct(x -> func.value(x).inner(e.getTangent().value(x)) * testFunctions.get(0).value(x.at(e.getTangentialDimension())), e);
		} else if (f != null)
		{
			return TPFaceIntegral.integrateNonTensorProduct(x ->
			{
				CoordinateVector val = func.value(x);
				if (f.getFlatDimension() == 0)
				{
					double testFunctionComponentValue =
						testFunctions.get(0).value(x.at(1)) * testFunctions.get(1).value(x.at(2));
					if (testComponent == 0)
						return val.at(2) * testFunctionComponentValue;
					else
						return -val.at(1) * testFunctionComponentValue;
				} else if (f.getFlatDimension() == 1)
				{
					double testFunctionComponentValue =
						testFunctions.get(0).value(x.at(0)) * testFunctions.get(1).value(x.at(2));
					if (testComponent == 0)
						return -val.at(2) * testFunctionComponentValue;
					else
						return val.at(0) * testFunctionComponentValue;
				} else
				{
					double testFunctionComponentValue =
						testFunctions.get(0).value(x.at(0)) * testFunctions.get(1).value(x.at(1));
					if (testComponent == 0)
						return val.at(1) * testFunctionComponentValue;
					else
						return -val.at(0) * testFunctionComponentValue;
				}
			}, f.getCell1Ds(), f.getFlatDimension(), f.getOtherCoordinate());
		} else
		{
			return TPCellIntegral.integrateNonTensorProduct(x ->
			{
				double testFunctionComponentValue =
					testFunctions.get(0).value(x.at(0)) * testFunctions.get(1).value(x.at(1)) * testFunctions.get(2).value(x.at(2));
				return func.value(x).at(testComponent) * testFunctionComponentValue;
			}, c.getCell1Ds());
		}
	}
	
	@Override
	public int compareTo(@NotNull NedelecNodeFuncional o)
	{
		if(getType() != o.getType())
			return Integer.compare(getType(), o.getType());
		else if(getType() == EDGETYPE)
			return getE().compareTo(o.getE());
		else if(getType() == FACETYPE)
			return getF().compareTo(o.getF());
		else
			return getC().compareTo(o.getC());
	}
}
