package basic;

import linalg.Tensor;
import linalg.Vector;

public abstract class CellIntegral<CT extends Cell<CT,?,?>, ST extends ShapeFunction<CT,
	?,?,?,?,?>>
{
	protected Function<?,?,?> weight;
	protected String name;
	
	protected CellIntegral()
	{
		weight = null;
		name = "NotReferencedByName";
	}
	public CellIntegral(Function<?,?,?> weight, String name)
	{
		this.weight = weight;
		this.name = name;
	}
	
	public CellIntegral(String name)
	{
		this.weight = ScalarFunction.constantFunction(1);
		this.name = name;
	}
	
	public abstract double evaluateCellIntegral(CT cell, ST shapeFunction1,
	                                            ST shapeFunction2);

}
