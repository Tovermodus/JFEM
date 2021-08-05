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


public class TaylorHoodStokesOrder
{
	public static void main(String[] args)
	{
		
		
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();CoordinateVector start = CoordinateVector.fromValues(0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		
		TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
				TPVectorCellIntegral.GRAD_GRAD);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(gradGrad);
		List<CellIntegral<TPCell,  QkQkFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		
		List<FaceIntegral< TPFace, QkQkFunction>> faceIntegrals = new ArrayList<>();
		
		MixedRightHandSideIntegral<TPCell,  ContinuousTPShapeFunction,
			ContinuousTPVectorFunction,QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(StokesReferenceSolution.rightHandSide(), TPVectorRightHandSideIntegral.VALUE));
		
		List<RightHandSideIntegral<TPCell,  QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		List<BoundaryRightHandSideIntegral< TPFace, QkQkFunction>> boundaryFaceIntegrals = new ArrayList<>();
		int polynomialDegree = 3;
		List<ScalarFunction> solutions = new ArrayList<>();
		List<VectorFunction> solutionsVec = new ArrayList<>();
		TaylorHoodSpace grid = null;
		for(int i = 0; i < 4; i++)
		{
			grid = new TaylorHoodSpace(start, end,
				Ints.asList(2*(int)Math.pow(2,i),2*(int)Math.pow(2,i)));
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
			grid.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
			IterativeSolver it = new IterativeSolver();
			Vector solution1 = it.solveCG(grid.getSystemMatrix(), grid.getRhs(), 1e-8);
			MixedFESpaceFunction<QkQkFunction> solut =
				new MixedFESpaceFunction<>(
					grid.getShapeFunctions(), solution1);
			solutions.add(solut.getPressureFunction());
			solutionsVec.add(solut.getVelocityFunction());
		}
		solutions.add(StokesReferenceSolution.pressureReferenceSolution());
		solutionsVec.add(StokesReferenceSolution.velocityReferenceSolution());
		System.out.println(ConvergenceOrderEstimator.estimateL20Scalar(solutions, grid.generatePlotPoints(50)));
		System.out.println(ConvergenceOrderEstimator.estimateL2Vector(solutionsVec,
			grid.generatePlotPoints(50)));
		
		
	}
}
