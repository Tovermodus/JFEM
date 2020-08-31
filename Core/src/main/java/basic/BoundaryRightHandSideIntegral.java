package basic;

public abstract class BoundaryRightHandSideIntegral<CT extends Cell<CT,FT>, FT extends Face<CT,FT>, ST extends ShapeFunction<CT,
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