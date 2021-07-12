package basic;

public abstract class FaceIntegral<FT extends Face<?,FT,?>, ST extends ShapeFunction<?,	FT,?,?,?,?>>
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