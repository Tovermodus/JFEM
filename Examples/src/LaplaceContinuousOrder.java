import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class LaplaceContinuousOrder
{
	public static void main(String[] args)
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		
		TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell,ContinuousTPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		ArrayList<FaceIntegral<TPFace,ContinuousTPVectorFunction>> faceIntegrals =
			new ArrayList<>();
		TPVectorRightHandSideIntegral<ContinuousTPVectorFunction> rightHandSideIntegral =
			new TPVectorRightHandSideIntegral<>(LaplaceReferenceSolution.vectorRightHandSide(),
				TPVectorRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell, ContinuousTPVectorFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		ArrayList<BoundaryRightHandSideIntegral<TPFace, ContinuousTPVectorFunction>> boundaryFaceIntegrals = new ArrayList<>();
		
		int polynomialDegree = 2;
		List<ScalarFunction> solutions = new ArrayList<>();
		List<VectorFunction> solutionsVec = new ArrayList<>();
		ContinuousTPFEVectorSpace grid = null;
		for(int i = 0; i < 5; i++)
		{
			grid = new ContinuousTPFEVectorSpace(start, end,
				Ints.asList(2*(int)Math.pow(2,i),2*(int)Math.pow(2,i)));
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
			grid.setBoundaryValues(LaplaceReferenceSolution.vectorBoundaryValues());
			IterativeSolver it = new IterativeSolver();
			Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-11);
			VectorFESpaceFunction<ContinuousTPVectorFunction> solut =
				new VectorFESpaceFunction<>(
					grid.getShapeFunctions(), solution1);
			solutionsVec.add(new VectorFunction()
			{
				@Override
				public int getRangeDimension()
				{
					return solut.getRangeDimension();
				}
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
