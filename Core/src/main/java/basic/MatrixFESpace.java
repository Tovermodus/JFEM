package basic;

import linalg.Matrix;
import linalg.Vector;

import java.util.List;

public interface MatrixFESpace<CT extends Cell<CT,FT>, FT extends  Face<CT,FT>,
	ST extends ShapeFunction<CT,FT,ST,valueT,gradientT,hessianT>,valueT,gradientT,hessianT,
	FST extends MatrixFESpace<CT,FT,ST,valueT,gradientT,hessianT,FST>> extends FESpace<CT,FT
	,ST,valueT,gradientT,hessianT,FST>, FESpaceTools<CT,FT,ST>
{
	void initializeSystemMatrix();
	
	void initializeRhs();
	
	Vector getRhs();
	
	Matrix getSystemMatrix();
	
	default void evaluateCellIntegrals(List<CellIntegral<CT, FT, ST>> cellIntegrals,
	                                   List<RightHandSideIntegral<CT, FT, ST>> rightHandSideIntegrals)
	{
		loopMatrixViaCell((K, u, v) ->
		{
			double integral = 0;
			for (CellIntegral<CT, FT, ST> cellIntegral :
				cellIntegrals)
			{
				integral += cellIntegral.evaluateCellIntegral(K, u, v);
			}
			return integral;
		}, this);
		loopRhsViaCell((K,  v) ->
		{
			double integral = 0;
			for (RightHandSideIntegral<CT, FT, ST> rightHandSideIntegral :
				rightHandSideIntegrals)
			{
				integral += rightHandSideIntegral.evaluateRightHandSideIntegral(K, v);
			}
			return integral;
		}, this);
	}
	
	default void evaluateFaceIntegrals(List<FaceIntegral<CT, FT, ST>> faceIntegrals,
	                                   List<BoundaryRightHandSideIntegral<CT, FT, ST>> boundaryRightHandSideIntegrals)
	{
		
		loopMatrixViaFace((F, u, v) ->
		{
			double integral = 0;
			for (FaceIntegral<CT, FT, ST> faceIntegral :
				faceIntegrals)
			{
				integral += faceIntegral.evaluateFaceIntegral(F, u, v);
			}
			return integral;
		}, this);
		loopRhsViaFace((F, v) ->
		{
			double integral = 0;
			for (BoundaryRightHandSideIntegral<CT, FT, ST> boundaryRightHandSideIntegral :
				boundaryRightHandSideIntegrals)
			{
				integral += boundaryRightHandSideIntegral.evaluateBoundaryRightHandSideIntegral(F, v);
			}
			return integral;
		}, this);
	}
}
