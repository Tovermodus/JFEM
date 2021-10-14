package basic;

import io.vavr.Function3;
import linalg.MutableMatrix;
import linalg.MutableVector;

import java.util.List;
import java.util.function.BiFunction;
import java.util.function.BiPredicate;

public interface MatrixFESpace<CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
	ST extends ShapeFunction<CT, FT, ?, ?, ?>> extends FESpace<CT, FT, ST>
{
	void initializeSystemMatrix();
	
	void initializeRhs();
	
	MutableVector getRhs();
	
	MutableMatrix getSystemMatrix();
	
	default void evaluateCellIntegrals(final List<CellIntegral<CT, ST>> cellIntegrals,
	                                   final List<RightHandSideIntegral<CT, ST>> rightHandSideIntegrals)
	{
		writeCellIntegralsToMatrix(cellIntegrals, getSystemMatrix());
		writeCellIntegralsToRhs(rightHandSideIntegrals, getRhs());
	}
	
	default void writeCellIntegralsToMatrix(final List<CellIntegral<CT, ST>> cellIntegrals,
	                                        final MutableMatrix s)
	{
		loopMatrixViaCell((K, u, v) ->
		                  {
			                  double integral = 0;
			                  for (final CellIntegral<CT, ST> cellIntegral :
				                  cellIntegrals)
			                  {
				                  integral += cellIntegral.evaluateCellIntegral(K, u, v);
			                  }
			                  return integral;
		                  }, s);
	}
	
	default void writeCellIntegralsToRhs(final List<RightHandSideIntegral<CT, ST>> rightHandSideIntegrals,
	                                     final MutableVector d)
	{
		loopRhsViaCell((K, v) ->
		               {
			               double integral = 0;
			               for (final RightHandSideIntegral<CT, ST> rightHandSideIntegral :
				               rightHandSideIntegrals)
			               {
				               integral += rightHandSideIntegral.evaluateRightHandSideIntegral(K, v);
			               }
			               return integral;
		               }, d);
	}
	
	default void writeCellIntegralsToRhs(final List<RightHandSideIntegral<CT, ST>> rightHandSideIntegrals,
	                                     final MutableVector d, final BiPredicate<CT, ST> functions)
	{
		loopRhsViaCell((K, v) ->
		               {
			               double integral = 0;
			               if (functions.test(K, v))
				               for (final RightHandSideIntegral<CT, ST> rightHandSideIntegral :
					               rightHandSideIntegrals)
				               {
					               integral += rightHandSideIntegral.evaluateRightHandSideIntegral(
						               K, v);
				               }
			               return integral;
		               }, d);
	}
	
	default void evaluateFaceIntegrals(final List<FaceIntegral<FT, ST>> faceIntegrals,
	                                   final List<BoundaryRightHandSideIntegral<FT, ST>> boundaryRightHandSideIntegrals)
	{
		writeFaceIntegralsToMatrix(faceIntegrals, getSystemMatrix());
		writeFaceIntegralsToRhs(boundaryRightHandSideIntegrals, getRhs());
	}
	
	default void writeFaceIntegralsToMatrix(final List<FaceIntegral<FT, ST>> faceIntegrals, final MutableMatrix s)
	{
		
		loopMatrixViaFace((F, u, v) ->
		                  {
			                  double integral = 0;
			                  for (final FaceIntegral<FT, ST> faceIntegral :
				                  faceIntegrals)
			                  {
				                  integral += faceIntegral.evaluateFaceIntegral(F, u, v);
			                  }
			                  return integral;
		                  }, s);
	}
	
	default void writeFaceIntegralsToRhs(final List<BoundaryRightHandSideIntegral<FT, ST>> boundaryRightHandSideIntegrals
		, final MutableVector d)
	{
		loopRhsViaFace((F, v) ->
		               {
			               double integral = 0;
			               for (final BoundaryRightHandSideIntegral<FT, ST> boundaryRightHandSideIntegral :
				               boundaryRightHandSideIntegrals)
			               {
				               integral +=
					               boundaryRightHandSideIntegral.evaluateBoundaryRightHandSideIntegral(
						               F, v);
			               }
			               return integral;
		               }, d);
	}
	
	default <T> void addToMatrix(final Function3<T, ST, ST, Double> integralEvaluation, final MutableMatrix s, final T K, final ST v,
	                             final ST u)
	{
		final double integral = integralEvaluation.apply(K, u, v);
		if (integral != 0)
			s.add(integral, v.getGlobalIndex(),
			      u.getGlobalIndex());
	}
	
	default <T> void addToVector(final BiFunction<T, ST, Double> integralEvaluation, final MutableVector d, final T K, final ST v)
	{
		final double integral = integralEvaluation.apply(K, v);
		if (integral != 0)
			d.add(integral, v.getGlobalIndex());
	}
	
	default void loopMatrixViaCell(final Function3<CT, ST, ST, Double> integralEvaluation, final MutableMatrix s)
	{
		forEachCell(K ->
		            {
			            forEachFunctionCombinationOnCell(K,
			                                             (u, v) -> addToMatrix(integralEvaluation, s, K, v,
			                                                                   u));
		            });
	}
	
	default void loopRhsViaCell(final BiFunction<CT, ST, Double> integralEvaluation, final MutableVector d)
	{
		forEachCell(K ->
			            forEachFunctionOnCell(K, (u) -> addToVector(integralEvaluation, d, K, u)));
	}
	
	default void loopMatrixViaFace(final Function3<FT, ST, ST, Double> integralEvaluation, final MutableMatrix s)
	{
		forEachFace(F ->
			            forEachFunctionCombinationOnFace(F,
			                                             (u, v) -> addToMatrix(integralEvaluation, s, F, v,
			                                                                   u)));
	}
	
	default void loopRhsViaFace(final BiFunction<FT, ST, Double> integralEvaluation, final MutableVector d)
	{
		forEachFace(F ->
			            forEachFunctionOnFace(F, (u) -> addToVector(integralEvaluation, d, F, u)));
	}
}
