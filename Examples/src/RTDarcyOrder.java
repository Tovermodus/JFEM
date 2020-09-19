import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.SparseMatrix;
import linalg.Vector;
import mixed.*;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;


public class RTDarcyOrder
{
	public static void main(String[] args)
	{
		
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		TPVectorCellIntegral<RTShapeFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		MixedCellIntegral<TPCell, TPFace,ContinuousTPShapeFunction, RTShapeFunction>
			divValue = new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, TPFace, ContinuousTPShapeFunction, RTShapeFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell, TPFace, MixedShapeFunction<TPCell, TPFace,ContinuousTPShapeFunction,
									RTShapeFunction>>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral<TPCell, TPFace, MixedShapeFunction<TPCell, TPFace,ContinuousTPShapeFunction,
					RTShapeFunction>>> faceIntegrals = new ArrayList<>();
		MixedRightHandSideIntegral<TPCell, TPFace, ContinuousTPShapeFunction, RTShapeFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromPressureIntegral(
				new TPRightHandSideIntegral<ContinuousTPShapeFunction>(
					LaplaceReferenceSolution.scalarRightHandSide(),TPRightHandSideIntegral.VALUE, false));
		List<RightHandSideIntegral<TPCell, TPFace, MixedShapeFunction<TPCell, TPFace,ContinuousTPShapeFunction,
					RTShapeFunction>>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		List<BoundaryRightHandSideIntegral<TPCell, TPFace, MixedShapeFunction<TPCell, TPFace,ContinuousTPShapeFunction,
					RTShapeFunction>>> boundaryFaceIntegrals = new ArrayList<>();
		MixedBoundaryRightHandSideIntegral<TPCell, TPFace, ContinuousTPShapeFunction, RTShapeFunction> dirichlet =
			MixedBoundaryRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorBoundaryFaceIntegral<RTShapeFunction>(LaplaceReferenceSolution.scalarBoundaryValues(),
					TPVectorBoundaryFaceIntegral.NORMAL_VALUE));
		boundaryFaceIntegrals.add(dirichlet);
		int polynomialDegree = 3;
		List<ScalarFunction> solutions = new ArrayList<>();
		List<VectorFunction> solutionsVec = new ArrayList<>();
		MixedRTSpace grid = null;
		for(int i = 0; i < 4; i++)
		{
			grid = new MixedRTSpace(start, end,
				Ints.asList(2*(int)Math.pow(2,i),2*(int)Math.pow(2,i)), polynomialDegree);
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
			IterativeSolver<SparseMatrix> it = new IterativeSolver<>();
			Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-9);
			MixedFESpaceFunction<TPCell,TPFace,ContinuousTPShapeFunction,RTShapeFunction> solut =
				new MixedFESpaceFunction<>(
					grid.getShapeFunctions(), solution1);
			solutions.add(new ScalarFunction()
			{
				@Override
				public int getDomainDimension()
				{
					return solut.getDomainDimension();
				}
				
				@Override
				public Double value(CoordinateVector pos)
				{
					return solut.value(pos).getPressure();
				}
			});
			solutionsVec.add(new VectorFunction()
			{
				@Override
				public int getDomainDimension()
				{
					return solut.getDomainDimension();
				}
				
				@Override
				public CoordinateVector value(CoordinateVector pos)
				{
					return solut.value(pos).getVelocity();
				}
			});
		}
		solutions.add(LaplaceReferenceSolution.scalarReferenceSolution());
		System.out.println(ConvergenceOrderEstimator.estimateL2Scalar(solutions, grid.generatePlotPoints(50)));
		
		
	}
}
