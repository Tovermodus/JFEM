package mixed;
import basic.*;
import com.google.common.collect.Lists;
import linalg.Matrix;
import linalg.Vector;

import java.util.List;
import java.util.function.DoubleFunction;
import java.util.function.DoubleSupplier;

public interface MixedFESpace<CT extends Cell<CT,FT>, FT extends Face<CT,FT>,ST extends MixedShapeFunction<CT,FT,ST,
	PF,VF>, PF extends ScalarShapeFunction<CT,FT,PF>,VF extends VectorShapeFunction<CT,FT,VF>> extends FESpace<CT
	, FT, ST, MixedValue, MixedGradient, MixedHessian, MixedFESpace<CT, FT, ST, PF, VF>>, FunctionSpaceTools<CT, FT, ST>
{
	void initializeSystemMatrix();
	void initializeRhs();
	
	Vector getRhs();
	
	Matrix getSystemMatrix();
	default void evaluatePressureCellIntegrals(){}
	default void evaluatePressureFaceIntegrals(){}
	default void evaluateVelocityCellIntegrals(){}
	default void evaluateVelocityFaceIntegrals(){}
	default void evaluateMixedCellIntegrals(){}
	default void evaluateMixedFaceIntegrals(){}
	
}
