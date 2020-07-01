package basic;

public abstract class FaceIntegral
{
	protected ScalarFunction weight;
	private TensorFunction vectorWeight;
	protected String name;
	public FaceIntegral(ScalarFunction weight, String name)
	{
		this.weight = weight;
		this.name = name;
	}
	public FaceIntegral(TensorFunction vectorWeight, String name)
	{
		this.vectorWeight = vectorWeight;
		this.name = name;
	}
	public FaceIntegral(String name)
	{
		this(ScalarFunction.oneFunction(), name);
	}
	
	public abstract double evaluateFaceIntegral(Face face, ScalarShapeFunction function1,
	                                     ScalarShapeFunction function2);

}
