import basic.CellIntegral;
import basic.FaceIntegral;
import basic.PerformanceArguments;
import basic.ScalarFunction;
import com.google.common.base.Stopwatch;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.jetbrains.annotations.NotNull;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class LaplaceContinuous
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		builder.build();
		
		System.out.println("output start");
		final CoordinateVector start = CoordinateVector.fromValues(-3, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 1;
		final double penalty = 100;
		final ContinuousTPFESpace grid2 = getContinuousTPFESpace(start, end, polynomialDegree, penalty);
		
		final CTPFESpace grid = getCtpfeSpace(start, end, polynomialDegree, penalty);
		System.out.println(grid.getSystemMatrix()
		                       .sub(grid2.getSystemMatrix())
		                       .absMaxElement());
//		final IterativeSolver it = new IterativeSolver();
//		final Vector solution1 = it.solveCG(grid.getSystemMatrix(), grid.getRhs(), 1e-3);
		System.out.println("solved");
//
//		final ScalarFESpaceFunction<CTPShapeFunction> solut =
//			new ScalarFESpaceFunction<>(
//				grid.getShapeFunctionMap(), solution1);
//		PlotWindow.addPlot(new ScalarPlot2D(solut, grid.generatePlotPoints(50), 50));
	}
	
	@NotNull
	private static ContinuousTPFESpace getContinuousTPFESpace(final CoordinateVector start,
	                                                          final CoordinateVector end,
	                                                          final int polynomialDegree,
	                                                          final double penalty)
	{
		final ContinuousTPFESpace grid2 = new ContinuousTPFESpace(start, end,
		                                                          new IntCoordinates(100, 100));
		final TPCellIntegral<ContinuousTPShapeFunction> gg2 =
			new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			                     TPCellIntegral.GRAD_GRAD);
		final ArrayList<CellIntegral<TPCell, ContinuousTPShapeFunction>> cellIntegrals2 =
			new ArrayList<>();
		cellIntegrals2.add(gg2);
		final ArrayList<FaceIntegral<TPFace, ContinuousTPShapeFunction>> faceIntegrals2 = new ArrayList<>();
		final TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral2 =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(4),
			                              TPRightHandSideIntegral.VALUE);
		grid2.assembleCells();
		grid2.assembleFunctions(polynomialDegree);
		grid2.initializeSystemMatrix();
		grid2.initializeRhs();
		Stopwatch s = Stopwatch.createStarted();
		System.out.println("integ start");
		grid2.evaluateCellIntegrals(cellIntegrals2, List.of(rightHandSideIntegral2));
		System.out.println(s.elapsed());
		s = Stopwatch.createStarted();
		grid2.evaluateFaceIntegrals(faceIntegrals2, List.of());
		System.out.println(s.elapsed());
		grid2.setBoundaryValues(ScalarFunction.constantFunction(0));
		return grid2;
	}
	
	@NotNull
	private static CTPFESpace getCtpfeSpace(final CoordinateVector start,
	                                        final CoordinateVector end,
	                                        final int polynomialDegree,
	                                        final double penalty)
	{
		Stopwatch s;
		final CTPFESpace grid = new CTPFESpace(start, end,
		                                       new IntCoordinates(100, 100));
		final TPCellIntegralOnReferenceCell gg =
			new TPCellIntegralOnReferenceCell(ScalarFunction.constantFunction(1),
			                                  TPCellIntegral.GRAD_GRAD, 2);
		
		final ArrayList<CellIntegral<TPCell, CTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		final ArrayList<FaceIntegral<TPFace, CTPShapeFunction>> faceIntegrals = new ArrayList<>();
		final TPRightHandSideIntegral<CTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(4),
			                              TPRightHandSideIntegral.VALUE);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		s = Stopwatch.createStarted();
		System.out.println("integ start");
		grid.evaluateCellIntegrals(cellIntegrals, List.of(rightHandSideIntegral));
		System.out.println(s.elapsed());
		s = Stopwatch.createStarted();
		grid.evaluateFaceIntegrals(faceIntegrals, List.of());
		System.out.println(s.elapsed());
		grid.setBoundaryValues(ScalarFunction.constantFunction(0));
		return grid;
	}
}
