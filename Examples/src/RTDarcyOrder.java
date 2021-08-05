import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;


public class RTDarcyOrder
{
	public static void main(String[] args)
	{
		
		
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		TPVectorCellIntegral<RTShapeFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction>
			divValue = new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell,  RTMixedFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral< TPFace, RTMixedFunction>> faceIntegrals = new ArrayList<>();
		MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromPressureIntegral(
				new TPRightHandSideIntegral<>(
					LaplaceReferenceSolution.scalarRightHandSide(),TPRightHandSideIntegral.VALUE));
		List<RightHandSideIntegral<TPCell,  RTMixedFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		List<BoundaryRightHandSideIntegral< TPFace, RTMixedFunction>> boundaryFaceIntegrals = new ArrayList<>();
		MixedBoundaryRightHandSideIntegral<TPFace, ContinuousTPShapeFunction, RTShapeFunction,
			RTMixedFunction> dirichlet =
			MixedBoundaryRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorBoundaryFaceIntegral<>(LaplaceReferenceSolution.scalarBoundaryValues(),
					TPVectorBoundaryFaceIntegral.NORMAL_VALUE));
		boundaryFaceIntegrals.add(dirichlet);
		int polynomialDegree = 3;
		List<ScalarFunction> solutions = new ArrayList<>();
		List<VectorFunction> solutionsVec = new ArrayList<>();
		MixedRTSpace grid = null;
		for(int i = 0; i < 3; i++)
		{
			grid = new MixedRTSpace(start, end,
				Ints.asList(2*(int)Math.pow(2,i),2*(int)Math.pow(2,i)));
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
			IterativeSolver it = new IterativeSolver();
			Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-9);
			MixedFESpaceFunction<RTMixedFunction> solut =
				new MixedFESpaceFunction<>(
					grid.getShapeFunctions(), solution1);
			solutions.add(solut.getPressureFunction());
			solutionsVec.add(solut.getVelocityFunction());
		}
		solutions.add(LaplaceReferenceSolution.scalarReferenceSolution());
		solutionsVec.add(LaplaceReferenceSolution.scalarReferenceSolution().getGradientFunction());
		System.out.println(ConvergenceOrderEstimator.estimateL2Scalar(solutions, grid.generatePlotPoints(50)));
		System.out.println(ConvergenceOrderEstimator.estimateL2Vector(solutionsVec,
			grid.generatePlotPoints(50)));
		
		
	}
}
