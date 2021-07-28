package basic;

import tensorproduct.QuadratureRule1D;

public abstract class RightHandSideIntegral<CT extends Cell<CT,?,?>,
	ST extends ShapeFunction<CT,
	?,?,?,?,?>>
{
	protected Function<?, ?, ?> rightHandSide;
	protected String name;
	public final QuadratureRule1D quadratureRule1D;
	public RightHandSideIntegral()
	{
		this.quadratureRule1D = QuadratureRule1D.Gauss5;
		rightHandSide = null;
		name = "NotReferencedByName";
	}
	public RightHandSideIntegral(Function<?, ?, ?> rightHandSide, String name)
	{
		this.rightHandSide = rightHandSide;
		this.name = name;
		this.quadratureRule1D = QuadratureRule1D.Gauss5;
	}
	public abstract double evaluateRightHandSideIntegral(CT cell, ST shapeFunction1);
	
}