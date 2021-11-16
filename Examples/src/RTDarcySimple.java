import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.MatrixPlot;
import linalg.Vector;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class RTDarcySimple
{
	public static void main(final String[] args)
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 1;
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
				new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(-1),
				                              TPRightHandSideIntegral.VALUE));
		final List<RightHandSideIntegral<TPCell, RTMixedFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		final List<BoundaryRightHandSideIntegral<TPFace, RTMixedFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		final MixedBoundaryRightHandSideIntegral<TPFace, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction>
			dirichlet =
			MixedBoundaryRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorBoundaryFaceIntegral<>(ScalarFunction.constantFunction(0),
				                                   TPVectorBoundaryFaceIntegral.NORMAL_VALUE));
		boundaryFaceIntegrals.add(dirichlet);
		
		final MixedRTSpace grid = new MixedRTSpace(start, end,
		                                           Ints.asList(3, 3));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		System.out.println("solve system: " + grid.getSystemMatrix()
		                                          .getRows() + "×" + grid.getSystemMatrix()
		                                                                 .getCols());
		
		final IterativeSolver i = new IterativeSolver();
		final Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-10);
		final MixedTPFESpaceFunction<RTMixedFunction> solut =
			new MixedTPFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		
		final PlotWindow p = new PlotWindow();
		p.addPlot(new MatrixPlot(grid.getSystemMatrix()));
		p.addPlot(new MixedPlot2D(solut, grid.generatePlotPoints(30), 30));
	}
}
