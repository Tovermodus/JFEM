package basic;

public abstract class BoundaryFaceIntegral
{
	protected ScalarFunction rightHandSide;
	protected ScalarFunction weight;
	public BoundaryFaceIntegral(ScalarFunction rightHandSide,ScalarFunction weight)
	{
		this.rightHandSide = rightHandSide;
		this.weight = weight;
	}

	public abstract double evaluateBoundaryFaceIntegral(Face f, ScalarShapeFunction v);
}
