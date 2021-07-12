package basic;

public interface NodeFunctional<valueT, gradientT, hessianT>
{
	FunctionSignature getFunctionSignature();
	double evaluate(Function<valueT, gradientT, hessianT> func);
	default boolean canEvaluate(FunctionSignature functionSignature)
	{
		 return functionSignature.getValueT().isAssignableFrom(functionSignature.getValueT())
			 && functionSignature.getGradientT().isAssignableFrom(functionSignature.getGradientT())
			 && functionSignature.getHessianT().isAssignableFrom(functionSignature.getHessianT());
	}
}
