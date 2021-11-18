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
import java.util.Map;

public class QkQkDarcy
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 1;
		final QkQkSpace grid = new QkQkSpace(start, end,
		                                     Ints.asList(10, 10));
		final TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue = new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
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
			MixedRightHandSideIntegral.fromPressureIntegral(
				new TPRightHandSideIntegral<>(
					LaplaceReferenceSolution.scalarRightHandSide(), TPRightHandSideIntegral.VALUE));
		final List<RightHandSideIntegral<TPCell, QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		final List<BoundaryRightHandSideIntegral<TPFace, QkQkFunction>> boundaryFaceIntegrals
			= new ArrayList<>();
		final MixedBoundaryRightHandSideIntegral<TPFace, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> dirichlet =
			MixedBoundaryRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorBoundaryFaceIntegral<>(LaplaceReferenceSolution.scalarBoundaryValues(),
				                                   TPVectorBoundaryFaceIntegral.NORMAL_VALUE));
		boundaryFaceIntegrals.add(dirichlet);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setPressureBoundaryValues(ScalarFunction.constantFunction(0));
		System.out.println("solve system: " + grid.getSystemMatrix()
		                                          .getRows() + "×" + grid.getSystemMatrix()
		                                                                 .getCols());
		//grid.A.makeParallelReady(12);
		if (grid.getRhs()
		        .getLength() < 100)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		final IterativeSolver i = new IterativeSolver();
		final Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-3);
		System.out.println("solved");
		System.out.println(solution1);
		System.out.println(grid.getSystemMatrix()
		                       .mvMul(solution1));
		final MixedTPFESpaceFunction<QkQkFunction> solut =
			new MixedTPFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		final ArrayList<Map<CoordinateVector, Double>> valList = new ArrayList<>();
		valList.add(solut.pressureValuesInPoints(grid.generatePlotPoints(50)));
		valList.add(solut.velocityComponentsInPoints(grid.generatePlotPoints(50), 0));
		valList.add(solut.velocityComponentsInPoints(grid.generatePlotPoints(50), 1));
		for (final ComposeMixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction, ContinuousTPVectorFunction> shapeFunction : grid.getShapeFunctions()
		                                                                                                                                .values())
		{
			if (shapeFunction.hasPressureFunction())
				valList.add(shapeFunction.pressureValuesInPoints(grid.generatePlotPoints(50)));
			if (shapeFunction.hasVelocityFunction())
			{
				valList.add(shapeFunction.velocityComponentsInPoints(grid.generatePlotPoints(50),
				                                                     shapeFunction.getVelocityFunction()
				                                                                  .getComponent()));
			}
		}
		new PlotFrame(valList, start, end);
	}
}
