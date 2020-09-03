package mixed;
import basic.*;
import linalg.Matrix;
import linalg.Vector;

import java.util.List;

public interface MixedFESpace<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,ST extends MixedShapeFunction<CT,FT,ST,
	PF,VF>, PF extends ScalarShapeFunction<CT,FT,PF>,VF extends VectorShapeFunction<CT,FT,VF>> extends MatrixFESpace<CT
	, FT, ST, MixedValue, MixedGradient, MixedHessian, MixedFESpace<CT, FT, ST, PF, VF>>, FESpaceTools<CT, FT,
	ST>
{	
}
