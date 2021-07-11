package basic;

public interface NodeFunctional<FT extends Function<valueT, gradientT, hessianT>, valueT, gradientT, hessianT>
{
	FunctionSignature getFunctionSignature();
	double evaluate(FT func);
	default boolean canEvaluate(FunctionSignature functionSignature)
	{
		 return functionSignature.getValueT().isAssignableFrom(functionSignature.getValueT())
			 && functionSignature.getGradientT().isAssignableFrom(functionSignature.getGradientT())
			 && functionSignature.getHessianT().isAssignableFrom(functionSignature.getHessianT());
	}
	
}
