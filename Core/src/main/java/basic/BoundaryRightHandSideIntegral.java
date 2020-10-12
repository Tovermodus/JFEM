package basic;

public abstract class BoundaryRightHandSideIntegral<FT extends Face<?,FT,?>, ST extends ShapeFunction<?,
	FT,?,ST,?,?,?>>
{
	protected Function<?, ?, ?> rightHandSide;
	protected String name;
	public BoundaryRightHandSideIntegral()
	{
		rightHandSide = null;
		name = "NotReferencedByName";
	}
	public BoundaryRightHandSideIntegral(Function<?, ?, ?> rightHandSide, String name)
	{
		this.rightHandSide = rightHandSide;
		this.name = name;
	}
	public abstract double evaluateBoundaryRightHandSideIntegral(FT face, ST shapeFunction1);
	
}