package basic;

public abstract class RightHandSideIntegral<CT extends Cell<CT,?,?>,
	ST extends ShapeFunction<CT,
	?,?,?,?,?,?>>
{
	protected Function<?, ?, ?> rightHandSide;
	protected String name;
	public RightHandSideIntegral()
	{
		rightHandSide = null;
		name = "NotReferencedByName";
	}
	public RightHandSideIntegral(Function<?, ?, ?> rightHandSide, String name)
	{
		this.rightHandSide = rightHandSide;
		this.name = name;
	}
	public abstract double evaluateRightHandSideIntegral(CT cell, ST shapeFunction1);
	
}