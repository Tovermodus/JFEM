package basic;

import tensorproduct.QuadratureRule1D;

public abstract class CellIntegral<CT extends Cell<CT,?>, ST extends ShapeFunction<CT,
	?, ?,?,?>>
{
	protected Function<?,?,?> weight;
	protected String name;
	
	public final QuadratureRule1D quadratureRule1D;
	protected CellIntegral(QuadratureRule1D quadratureRule1D)
	{
		this.quadratureRule1D = quadratureRule1D;
		weight = null;
		name = "NotReferencedByName";
	}
	public CellIntegral(Function<?, ?, ?> weight, String name, QuadratureRule1D quadratureRule1D)
	{
		this.weight = weight;
		this.name = name;
		this.quadratureRule1D = quadratureRule1D;
	}
	
	public CellIntegral(String name, QuadratureRule1D quadratureRule1D)
	{
		this.quadratureRule1D = quadratureRule1D;
		this.weight = ScalarFunction.constantFunction(1);
		this.name = name;
	}
	
	public abstract double evaluateCellIntegral(CT cell, ST shapeFunction1,
	                                            ST shapeFunction2);

}
