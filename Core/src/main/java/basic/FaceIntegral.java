package basic;

import tensorproduct.QuadratureRule1D;

public abstract class FaceIntegral<FT extends Face<?,FT,?>, ST extends ShapeFunction<?,	FT,?,?,?,?>>
{
	protected Function<?, ?, ?> weight;
	protected String name;
	public final QuadratureRule1D quadratureRule1D;
	public FaceIntegral(QuadratureRule1D quadratureRule1D)
	{
		this.quadratureRule1D = quadratureRule1D;
		weight = null;
		name = "NotReferencedByName";
	}
	public FaceIntegral(Function<?, ?, ?> weight, String name, QuadratureRule1D quadratureRule1D)
	{
		this.weight = weight;
		this.name = name;
		this.quadratureRule1D = quadratureRule1D;
	}
	public FaceIntegral(String name, QuadratureRule1D quadratureRule1D)
	{
		this.quadratureRule1D = quadratureRule1D;
		this.weight = ScalarFunction.constantFunction(1);
		this.name = name;
	}
	public abstract double evaluateFaceIntegral(FT face, ST shapeFunction1,
	                                            ST shapeFunction2);
	
}