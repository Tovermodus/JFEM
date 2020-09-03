package basic;

public abstract class FaceIntegral<CT extends Cell<CT,FT>, FT extends Face<CT,FT>, ST extends ShapeFunction<CT,
	FT,ST,?,?,?>>
{
	protected Function<?, ?, ?> weight;
	protected String name;
	public FaceIntegral()
	{
		weight = null;
		name = "NotReferencedByName";
	}
	public FaceIntegral(Function<?, ?, ?> weight, String name)
	{
		this.weight = weight;
		this.name = name;
	}
	public FaceIntegral(String name)
	{
		this.weight = ScalarFunction.constantFunction(1);
		this.name = name;
	}
	public abstract double evaluateFaceIntegral(FT face, ST shapeFunction1,
	                                            ST shapeFunction2);
	
}