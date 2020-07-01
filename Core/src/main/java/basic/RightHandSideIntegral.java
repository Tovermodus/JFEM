package basic;

public abstract class RightHandSideIntegral
{
	protected ScalarFunction rightHandSide;
	protected ScalarFunction weight;
	public RightHandSideIntegral(ScalarFunction rightHandSide, ScalarFunction weight)
	{
		this.rightHandSide = rightHandSide;
		this.weight = weight;
	}

	public RightHandSideIntegral(ScalarFunction rightHandSide)
	{
		this(rightHandSide, ScalarFunction.oneFunction());
	}

	public abstract double evaluateRightHandSideIntegral(Cell k, ScalarShapeFunction v);
	
}
