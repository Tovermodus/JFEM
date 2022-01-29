import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import mixed.*;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class QkQkMaxwellOrder
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(0, 0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1, 1);
		final int polynomialDegree = 2;
		final TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.ROT_ROT);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vv
			=
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		final List<CellIntegral<TPCell, QkQkFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		final List<FaceIntegral<TPFace, QkQkFunction>> faceIntegrals = new ArrayList<>();
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(MaxwellReferenceSolution.rightHandSide(),
				                                    TPVectorRightHandSideIntegral.VALUE));
		final List<RightHandSideIntegral<TPCell, QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		final List<BoundaryRightHandSideIntegral<TPFace, QkQkFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		QkQkSpace grid = null;
		final List<ScalarFunction> solutions = new ArrayList<>();
		final List<VectorFunction> solutionsVec = new ArrayList<>();
		for (int i = 0; i < 3; i++)
		{
			grid = new QkQkSpace(start, end,
			                     Ints.asList(3 * (int) Math.pow(2, i),
			                                 3 * (int) Math.pow(2, i),
			                                 3 * (int) Math.pow(2, i)));
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			System.out.println("velocity Boundary");
			grid.setVelocityBoundaryValues(MaxwellReferenceSolution.vectorBoundaryValues());
			System.out.println("pressure Boundary");
			grid.setPressureBoundaryValues(MaxwellReferenceSolution.pressureBoundaryValues());
			System.out.println("Cell Integrals");
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			System.out.println("Face Integrals");
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
			System.out.println("solve system: " + grid.getSystemMatrix()
			                                          .getRows() + "×" + grid.getSystemMatrix()
			                                                                 .getCols());
			final IterativeSolver it = new IterativeSolver();
			final Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-7);
			final MixedTPFESpaceFunction<QkQkFunction> solut =
				new MixedTPFESpaceFunction<>(
					grid.getShapeFunctionMap(), solution1);
			solutions.add(solut.getPressureFunction());
			solutionsVec.add(solut.getVelocityFunction());
		}
		solutions.add(MaxwellReferenceSolution.pressureReferenceSolution());
		solutionsVec.add(MaxwellReferenceSolution.velocityReferenceSolution());
		System.out.println(ConvergenceOrderEstimator.estimateL2Vector(solutionsVec,
		                                                              grid.generatePlotPoints(20)));
		System.out.println(ConvergenceOrderEstimator.estimateL2Scalar(solutions, grid.generatePlotPoints(20)));
	}
}
