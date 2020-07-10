package basic;

public abstract class RightHandSideIntegral<CT extends Cell<CT,FT,ST>, FT extends Face<CT,FT,ST>, ST extends ShapeFunction<CT,
	FT,ST,?,?,?>>
{
	protected Function<?, ?, ?> rightHandSide;
	protected String name;
	public RightHandSideIntegral(Function<?, ?, ?> rightHandSide, String name)
	{
		this.rightHandSide = rightHandSide;
		this.name = name;
	}
	public abstract double evaluateRightHandSideIntegral(CT cell, ST shapeFunction1);
	
}