package basic;


import java.util.Set;

import linalg.*;
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
	
	MutableVector getRhs();
	
	MutableMatrix getSystemMatrix();
	default Set<Integer> getFixedNodeIndices()
	{
		return new TreeSet<Integer>();
	}
	default void evaluateCellIntegrals(List<CellIntegral<CT, ST>> cellIntegrals,
	                                   List<RightHandSideIntegral<CT, ST>> rightHandSideIntegrals)
	{
		writeCellIntegralsToMatrix(cellIntegrals, getSystemMatrix());
		writeCellIntegralsToRhs(rightHandSideIntegrals, getRhs());
	}
	default void writeCellIntegralsToMatrix(List<CellIntegral<CT, ST>> cellIntegrals,
	                                   MutableMatrix s)
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
		}, this, s);
	}
	default void writeCellIntegralsToRhs(List<RightHandSideIntegral<CT, ST>> rightHandSideIntegrals,
	                                   MutableVector d)
	{
		loopRhsViaCell((K,  v) ->
		{
			double integral = 0;
			for (RightHandSideIntegral<CT, ST> rightHandSideIntegral :
				rightHandSideIntegrals)
			{
				integral += rightHandSideIntegral.evaluateRightHandSideIntegral(K, v);
			}
			return integral;
		}, this, d);
	}
	
	default void evaluateFaceIntegrals(List<FaceIntegral<FT, ST>> faceIntegrals,
	                                   List<BoundaryRightHandSideIntegral<FT, ST>> boundaryRightHandSideIntegrals)
	{
		writeFaceIntegralsToMatrix(faceIntegrals, getSystemMatrix());
		writeFaceIntegralsToRhs(boundaryRightHandSideIntegrals, getRhs());
	}
	default void writeFaceIntegralsToMatrix(List<FaceIntegral<FT, ST>> faceIntegrals, MutableMatrix s)
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
		}, this, s);
	}
	default void writeFaceIntegralsToRhs(List<BoundaryRightHandSideIntegral<FT, ST>> boundaryRightHandSideIntegrals
		, MutableVector d)
	{
		loopRhsViaFace((F, v) ->
		{
			double integral = 0;
			for (BoundaryRightHandSideIntegral<FT, ST> boundaryRightHandSideIntegral :
				boundaryRightHandSideIntegrals)
			{
				integral += boundaryRightHandSideIntegral.evaluateBoundaryRightHandSideIntegral(F, v);
			}
			return integral;
		}, this, d);
	}
}
