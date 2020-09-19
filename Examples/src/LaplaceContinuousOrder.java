import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.SparseMatrix;
import linalg.Vector;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

public class LaplaceContinuousOrder
{
	public static void main(String[] args)
	{
		
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		
		TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		ArrayList<FaceIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> faceIntegrals = new ArrayList<>();
		TPVectorRightHandSideIntegral<ContinuousTPVectorFunction> rightHandSideIntegral =
			new TPVectorRightHandSideIntegral<>(LaplaceReferenceSolution.vectorRightHandSide(),
				TPVectorRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		ArrayList<BoundaryRightHandSideIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> boundaryFaceIntegrals = new ArrayList<>();
		
		int polynomialDegree = 2;
		List<ScalarFunction> solutions = new ArrayList<>();
		List<VectorFunction> solutionsVec = new ArrayList<>();
		ContinuousTPFEVectorSpace grid = null;
		for(int i = 0; i < 5; i++)
		{
			grid = new ContinuousTPFEVectorSpace(start, end,
				Ints.asList(2*(int)Math.pow(2,i),2*(int)Math.pow(2,i)), polynomialDegree);
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
			grid.setBoundaryValues(LaplaceReferenceSolution.vectorBoundaryValues());
			IterativeSolver<SparseMatrix> it = new IterativeSolver<>();
			Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-11);
			VectorFESpaceFunction<ContinuousTPVectorFunction> solut =
				new VectorFESpaceFunction<>(
					grid.getShapeFunctions(), solution1);
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
					return solut.value(pos);
				}
			});
		}
		solutionsVec.add(LaplaceReferenceSolution.vectorReferenceSolution());
		System.out.println(ConvergenceOrderEstimator.estimateL2Vector(solutionsVec,
			grid.generatePlotPoints(40)));
		
	}
}
