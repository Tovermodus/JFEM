package basic;

public abstract class FaceIntegral<CT extends Cell<CT,FT,ST>, FT extends Face<CT,FT,ST>, ST extends ShapeFunction<CT,
	FT,ST,?,?,?>>
{
	protected Function<?, ?, ?> weight;
	protected String name;
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