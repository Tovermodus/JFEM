import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class RTDarcyOrder
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final TPVectorCellIntegral<RTShapeFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction>
			divValue = new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		final List<CellIntegral<TPCell, RTMixedFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		final List<FaceIntegral<TPFace, RTMixedFunction>> faceIntegrals = new ArrayList<>();
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction>
			rightHandSideIntegral =
			MixedRightHandSideIntegral.fromPressureIntegral(
				new TPRightHandSideIntegral<>(
					LaplaceReferenceSolution.scalarRightHandSide(), TPRightHandSideIntegral.VALUE));
		final List<RightHandSideIntegral<TPCell, RTMixedFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		final List<BoundaryRightHandSideIntegral<TPFace, RTMixedFunction>> boundaryFaceIntegrals = new ArrayList<>();
		final MixedBoundaryRightHandSideIntegral<TPFace, ContinuousTPShapeFunction, RTShapeFunction,
			RTMixedFunction> dirichlet =
			MixedBoundaryRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorBoundaryFaceIntegral<>(LaplaceReferenceSolution.scalarBoundaryValues(),
				                                   TPVectorBoundaryFaceIntegral.NORMAL_VALUE));
		boundaryFaceIntegrals.add(dirichlet);
		final int polynomialDegree = 3;
		final List<ScalarFunction> solutions = new ArrayList<>();
		final List<VectorFunction> solutionsVec = new ArrayList<>();
		MixedRTSpace grid = null;
		for (int i = 0; i < 3; i++)
		{
			grid = new MixedRTSpace(start, end,
			                        Ints.asList(2 * (int) Math.pow(2, i), 2 * (int) Math.pow(2, i)));
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
			final IterativeSolver it = new IterativeSolver();
			final Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-9);
			final MixedTPFESpaceFunction<RTMixedFunction> solut =
				new MixedTPFESpaceFunction<>(
					grid.getShapeFunctions(), solution1);
			solutions.add(solut.getPressureFunction());
			solutionsVec.add(solut.getVelocityFunction());
		}
		solutions.add(LaplaceReferenceSolution.scalarReferenceSolution());
		solutionsVec.add(LaplaceReferenceSolution.scalarReferenceSolution()
		                                         .getGradientFunction());
		System.out.println(ConvergenceOrderEstimator.estimateL2Scalar(solutions, grid.generatePlotPoints(50)));
		System.out.println(ConvergenceOrderEstimator.estimateL2Vector(solutionsVec,
		                                                              grid.generatePlotPoints(50)));
	}
}
