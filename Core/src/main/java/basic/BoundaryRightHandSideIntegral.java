package basic;

public abstract class BoundaryRightHandSideIntegral<CT extends Cell<CT,FT,ST>, FT extends Face<CT,FT,ST>, ST extends ShapeFunction<CT,
	FT,ST,?,?,?>>
{
	protected Function<?, ?, ?> rightHandSide;
	protected String name;
	public BoundaryRightHandSideIntegral(Function<?, ?, ?> rightHandSide, String name)
	{
		this.rightHandSide = rightHandSide;
		this.name = name;
	}
	public abstract double evaluateBoundaryRightHandSideIntegral(FT face, ST shapeFunction1);
	
}