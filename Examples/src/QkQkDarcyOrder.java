import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.*;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;


public class QkQkDarcyOrder
{
	public static void main(String[] args)
	{
		
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		MixedCellIntegral<TPCell, TPFace,ContinuousTPShapeFunction, ContinuousTPVectorFunction>
			divValue = new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, TPFace, ContinuousTPShapeFunction, ContinuousTPVectorFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell, TPFace, MixedShapeFunction<TPCell, TPFace,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral<TPCell, TPFace, MixedShapeFunction<TPCell, TPFace,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> faceIntegrals = new ArrayList<>();
		MixedRightHandSideIntegral<TPCell, TPFace, ContinuousTPShapeFunction, ContinuousTPVectorFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromPressureIntegral(
				new TPRightHandSideIntegral<ContinuousTPShapeFunction>(
					LaplaceReferenceSolution.scalarRightHandSide(),TPRightHandSideIntegral.VALUE, false));
		List<RightHandSideIntegral<TPCell, TPFace, MixedShapeFunction<TPCell, TPFace,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		List<BoundaryRightHandSideIntegral<TPCell, TPFace, MixedShapeFunction<TPCell, TPFace,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> boundaryFaceIntegrals = new ArrayList<>();
		MixedBoundaryRightHandSideIntegral<TPCell, TPFace, ContinuousTPShapeFunction, ContinuousTPVectorFunction> dirichlet =
			MixedBoundaryRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorBoundaryFaceIntegral<ContinuousTPVectorFunction>(LaplaceReferenceSolution.scalarBoundaryValues(),
					TPVectorBoundaryFaceIntegral.NORMAL_VALUE));
		boundaryFaceIntegrals.add(dirichlet);
		int polynomialDegree = 1;
		List<ScalarFunction> solutions = new ArrayList<>();
		List<VectorFunction> solutionsVec = new ArrayList<>();
		QkQkSpace grid = null;
		for(int i = 0; i < 3; i++)
		{
			grid = new QkQkSpace(start, end,
				Ints.asList(2*(int)Math.pow(2,i),2*(int)Math.pow(2,i)), polynomialDegree);
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
			IterativeSolver<SparseMatrix> it = new IterativeSolver<>();
			Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-9);
			MixedFESpaceFunction<TPCell,TPFace,ContinuousTPShapeFunction,ContinuousTPVectorFunction> solut =
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
