package basic;


import io.vavr.Function3;
import linalg.*;

import java.util.List;
import java.util.function.BiFunction;

public interface MatrixFESpace<CT extends Cell<CT,FT>, FT extends  Face<CT,FT>,
	ST extends ShapeFunction<CT,FT, ?,?,?>> extends FESpace<CT,FT, ST>
{
	void initializeSystemMatrix();
	
	void initializeRhs();
	
	MutableVector getRhs();
	
	MutableMatrix getSystemMatrix();
	
	default void evaluateCellIntegrals(List<CellIntegral<CT,ST>> cellIntegrals,
	                                   List<RightHandSideIntegral<CT,ST>> rightHandSideIntegrals)
	{
		writeCellIntegralsToMatrix(cellIntegrals, getSystemMatrix());
		writeCellIntegralsToRhs(rightHandSideIntegrals, getRhs());
	}
	default void writeCellIntegralsToMatrix(List<CellIntegral<CT,ST>> cellIntegrals,
	                                   MutableMatrix s)
	{
		loopMatrixViaCell((K, u, v) ->
		{
			double integral = 0;
			for (CellIntegral<CT,ST> cellIntegral :
				cellIntegrals)
			{
				integral += cellIntegral.evaluateCellIntegral(K, u, v);
			}
			return integral;
		},  s);
	}
	default void writeCellIntegralsToRhs(List<RightHandSideIntegral<CT,ST>> rightHandSideIntegrals,
	                                   MutableVector d)
	{
		loopRhsViaCell((K,  v) ->
		{
			double integral = 0;
			for (RightHandSideIntegral<CT,ST> rightHandSideIntegral :
				rightHandSideIntegrals)
			{
				integral += rightHandSideIntegral.evaluateRightHandSideIntegral(K, v);
			}
			return integral;
		},  d);
	}
	
	default void evaluateFaceIntegrals(List<FaceIntegral<FT,ST>> faceIntegrals,
	                                   List<BoundaryRightHandSideIntegral<FT,ST>> boundaryRightHandSideIntegrals)
	{
		writeFaceIntegralsToMatrix(faceIntegrals, getSystemMatrix());
		writeFaceIntegralsToRhs(boundaryRightHandSideIntegrals, getRhs());
	}
	default void writeFaceIntegralsToMatrix(List<FaceIntegral<FT,ST>> faceIntegrals, MutableMatrix s)
	{
		
		loopMatrixViaFace((F, u, v) ->
		{
			double integral = 0;
			for (FaceIntegral<FT,ST> faceIntegral :
				faceIntegrals)
			{
				integral += faceIntegral.evaluateFaceIntegral(F, u, v);
			}
			return integral;
		}, s);
	}
	default void writeFaceIntegralsToRhs(List<BoundaryRightHandSideIntegral<FT,ST>> boundaryRightHandSideIntegrals
		, MutableVector d)
	{
		loopRhsViaFace((F, v) ->
		{
			double integral = 0;
			for (BoundaryRightHandSideIntegral<FT,ST> boundaryRightHandSideIntegral :
				boundaryRightHandSideIntegrals)
			{
				integral += boundaryRightHandSideIntegral.evaluateBoundaryRightHandSideIntegral(F, v);
			}
			return integral;
		},  d);
	}
	
	default <T> void addToMatrix(Function3<T, ST, ST, Double> integralEvaluation, MutableMatrix s, T K, ST v,
	                             ST u)
	{
		double integral = integralEvaluation.apply(K, u, v);
		if (integral != 0)
			s.add(integral, v.getGlobalIndex(),
				u.getGlobalIndex());
	}
	default <T> void addToVector(BiFunction<T, ST, Double> integralEvaluation, MutableVector d, T K, ST v)
	{
		double integral = integralEvaluation.apply(K, v);
		if (integral != 0)
			d.add(integral, v.getGlobalIndex());
	}
	
	default void loopMatrixViaCell(Function3<CT, ST, ST, Double> integralEvaluation, MutableMatrix s)
	{
		forEachCell(K ->
			forEachFunctionCombinationOnCell(K, (u,v)-> addToMatrix(integralEvaluation, s, K, v, u)));
	}
	
	default void loopRhsViaCell(BiFunction<CT, ST, Double> integralEvaluation,MutableVector d)
	{
		forEachCell(K ->
			forEachFunctionOnCell(K, (u)-> addToVector(integralEvaluation, d, K, u)));
	}
	
	
	default void loopMatrixViaFace(Function3<FT, ST, ST, Double> integralEvaluation, MutableMatrix s)
	{
		forEachFace(F ->
			forEachFunctionCombinationOnFace(F, (u,v)-> addToMatrix(integralEvaluation, s, F, v, u)));
	}
	
	default void loopRhsViaFace(BiFunction<FT, ST, Double> integralEvaluation, MutableVector d)
	{
		forEachFace(F ->
			forEachFunctionOnFace(F, (u)-> addToVector(integralEvaluation, d, F, u)));
	}
	
	
}
