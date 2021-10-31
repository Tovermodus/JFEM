import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import tensorproduct.ContinuousTPFEVectorSpace;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class LaplaceContinuousOrder
{
	public static void main(final String[] args)
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		
		final TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
		final ArrayList<CellIntegral<TPCell, ContinuousTPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		final ArrayList<FaceIntegral<TPFace, ContinuousTPVectorFunction>> faceIntegrals =
			new ArrayList<>();
		final TPVectorRightHandSideIntegral<ContinuousTPVectorFunction> rightHandSideIntegral =
			new TPVectorRightHandSideIntegral<>(LaplaceReferenceSolution.vectorRightHandSide(),
			                                    TPVectorRightHandSideIntegral.VALUE);
		final ArrayList<RightHandSideIntegral<TPCell, ContinuousTPVectorFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		final ArrayList<BoundaryRightHandSideIntegral<TPFace, ContinuousTPVectorFunction>> boundaryFaceIntegrals
			= new ArrayList<>();
		
		final int polynomialDegree = 2;
		final List<ScalarFunction> solutions = new ArrayList<>();
		final List<VectorFunction> solutionsVec = new ArrayList<>();
		ContinuousTPFEVectorSpace grid = null;
		for (int i = 0; i < 4; i++)
		{
			grid = new ContinuousTPFEVectorSpace(start, end,
			                                     Ints.asList(2 * (int) Math.pow(2, i),
			                                                 2 * (int) Math.pow(2, i)));
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
			grid.setBoundaryValues(LaplaceReferenceSolution.vectorBoundaryValues());
			final IterativeSolver it = new IterativeSolver();
			final Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-11);
			final VectorFESpaceFunction<ContinuousTPVectorFunction> solut =
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
				public CoordinateVector value(final CoordinateVector pos)
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
