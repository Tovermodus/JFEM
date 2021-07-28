package basic;

import tensorproduct.QuadratureRule1D;

public abstract class BoundaryRightHandSideIntegral<FT extends Face<?,FT,?>, ST extends ShapeFunction<?,
	FT,?,?,?,?>>
{
	protected Function<?, ?, ?> rightHandSide;
	protected String name;
	public final QuadratureRule1D quadratureRule1D;
	public BoundaryRightHandSideIntegral()
	{
		this.quadratureRule1D = QuadratureRule1D.Gauss5;
		rightHandSide = null;
		name = "NotReferencedByName";
	}
	public BoundaryRightHandSideIntegral(Function<?, ?, ?> rightHandSide, String name)
	{
		this.rightHandSide = rightHandSide;
		this.name = name;
		this.quadratureRule1D = QuadratureRule1D.Gauss5;
	}
	public abstract double evaluateBoundaryRightHandSideIntegral(FT face, ST shapeFunction1);
	
}