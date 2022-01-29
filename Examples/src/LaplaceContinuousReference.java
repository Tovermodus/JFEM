import basic.PerformanceArguments;
import basic.PlotWindow;
import basic.ScalarFESpaceFunction;
import basic.ScalarPlot2D;
import com.google.common.primitives.Ints;
import linalg.*;
import tensorproduct.ContinuousTPFESpace;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.TPCellIntegral;
import tensorproduct.TPRightHandSideIntegral;

import java.util.List;

public class LaplaceContinuousReference
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		
		final TPCellIntegral<ContinuousTPShapeFunction> gg =
			new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		final TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
			                              TPRightHandSideIntegral.VALUE);
		final int polynomialDegree = 1;
		final ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end,
		                                                         Ints.asList(3 * 64, 3 * 64));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		final SparseMatrix s = new SparseMatrix(grid.getShapeFunctionMap()
		                                            .size(),
		                                        grid.getShapeFunctionMap()
		                                            .size());
		final DenseVector d = new DenseVector(grid.getShapeFunctionMap()
		                                          .size());
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.writeCellIntegralsToMatrix(List.of(gg), s);
		grid.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), d);
		grid.writeBoundaryValuesTo(LaplaceReferenceSolution.scalarBoundaryValues(), s, d);
		final IterativeSolver it = new IterativeSolver();
		final Vector solution1 = it.solveGMRES(s, d, 1e-11);
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> solut =
			new ScalarFESpaceFunction<>(
				grid.getShapeFunctionMap(), solution1);
		PlotWindow.addPlot(new ScalarPlot2D(solut, grid.generatePlotPoints(30), 30));
		PlotWindow.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.scalarReferenceSolution(),
		                                    grid.generatePlotPoints(30),
		                                    30));
	}
}
