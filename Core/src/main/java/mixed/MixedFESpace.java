package mixed;
import basic.*;

public interface MixedFESpace<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,
	PF extends ScalarShapeFunction<CT
	,FT>,VF extends VectorShapeFunction<CT,FT>> extends MatrixFESpace<CT
	, FT, MixedShapeFunction<CT,FT, PF,VF>>, FESpaceTools<CT, FT,
	MixedShapeFunction<CT,FT, PF,VF>>, Assembleable
{

}
