package basic;


import java.util.Set;
import linalg.Matrix;
import linalg.Vector;
import mixed.MixedCellIntegral;
import mixed.MixedTPCellIntegral;

import java.util.List;
import java.util.TreeSet;

public interface MatrixFESpace<CT extends Cell<CT,FT,ET>, FT extends  Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>,
	ST extends ShapeFunction<CT,FT,ET,ST,valueT,gradientT,hessianT>,valueT,gradientT,hessianT
	> extends FESpace<CT,FT,ET
	,ST,valueT,gradientT,hessianT, MatrixFESpace<CT,FT,ET,ST,valueT,gradientT,hessianT>>, FESpaceTools<CT,FT,ET,ST>
{
	void initializeSystemMatrix();
	
	void initializeRhs();
	
	Vector getRhs();
	
	Matrix getSystemMatrix();
	default Set<Integer> getFixedNodeIndices()
	{
		return new TreeSet<Integer>();
	}
	default void evaluateCellIntegrals(List<CellIntegral<CT, ST>> cellIntegrals,
	                                   List<RightHandSideIntegral<CT, ST>> rightHandSideIntegrals)
	{
		loopMatrixViaCell((K, u, v) ->
		{
			double integral = 0;
			for (CellIntegral<CT, ST> cellIntegral :
				cellIntegrals)
			{
				integral += cellIntegral.evaluateCellIntegral(K, u, v);
			}
			return integral;
		}, this);
		loopRhsViaCell((K,  v) ->
		{
			double integral = 0;
			for (RightHandSideIntegral<CT, ST> rightHandSideIntegral :
				rightHandSideIntegrals)
			{
				integral += rightHandSideIntegral.evaluateRightHandSideIntegral(K, v);
			}
			return integral;
		}, this);
	}
	
	default void evaluateFaceIntegrals(List<FaceIntegral<FT, ST>> faceIntegrals,
	                                   List<BoundaryRightHandSideIntegral<FT, ST>> boundaryRightHandSideIntegrals)
	{
		
		loopMatrixViaFace((F, u, v) ->
		{
			double integral = 0;
			for (FaceIntegral<FT, ST> faceIntegral :
				faceIntegrals)
			{
				integral += faceIntegral.evaluateFaceIntegral(F, u, v);
			}
			return integral;
		}, this);
		loopRhsViaFace((F, v) ->
		{
			double integral = 0;
			for (BoundaryRightHandSideIntegral<FT, ST> boundaryRightHandSideIntegral :
				boundaryRightHandSideIntegrals)
			{
				integral += boundaryRightHandSideIntegral.evaluateBoundaryRightHandSideIntegral(F, v);
			}
			return integral;
		}, this);
	}
}
