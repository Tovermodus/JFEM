package basic;

public interface NodeFunctional<FT extends Function<valueT, gradientT, hessianT>, valueT, gradientT, hessianT>
{
	
	double evaluate(FT func);
}
