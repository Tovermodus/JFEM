package basic;

public interface NodeFunctionalShapeFunction<CT extends Cell<CT,FT, ET>, FT extends Face<CT,FT, ET>,
	ET extends Edge<CT,FT,ET>, ST extends ShapeFunction<CT,FT,ET, ST,
	valueT,gradientT,hessianT>, valueT,	gradientT,
	hessianT> extends ShapeFunction<CT,FT,ET,ST,valueT,gradientT,hessianT>
{

}
