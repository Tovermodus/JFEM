package mixed;
import basic.*;
import linalg.Matrix;
import linalg.Vector;

import java.util.List;

public interface MixedFESpace<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>,
	ET extends Edge<CT,FT,ET>,PF extends ScalarShapeFunction<CT
	,FT,ET,
	PF>,VF extends VectorShapeFunction<CT,FT,ET,VF>> extends MatrixFESpace<CT
	, FT, ET,MixedShapeFunction<CT,FT,ET,PF,VF>, MixedValue, MixedGradient, MixedHessian>, FESpaceTools<CT, FT,ET,
	MixedShapeFunction<CT,FT,ET,PF,VF>>, Assembleable
{

}
