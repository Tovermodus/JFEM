package mixed;
import basic.*;
import linalg.Matrix;
import linalg.Vector;

import java.util.List;

public interface MixedFESpace<CT extends Cell<CT,FT>, FT extends Face<CT,FT>, PF extends ScalarShapeFunction<CT,FT,PF>,VF extends VectorShapeFunction<CT,FT,VF>> extends MatrixFESpace<CT
	, FT, MixedShapeFunction<CT,FT,PF,VF>, MixedValue, MixedGradient, MixedHessian>, FESpaceTools<CT, FT,
	MixedShapeFunction<CT,FT,PF,VF>>
{
}
