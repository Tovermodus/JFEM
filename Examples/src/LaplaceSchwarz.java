import basic.*;
import linalg.*;
import schwarz.AdditiveSubspaceCorrection;
import schwarz.CGSolver;
import schwarz.CartesianUpFrontSchwarz;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;

import java.util.ArrayList;

public class LaplaceSchwarz
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		
		final long startTime = System.nanoTime();
		
		System.out.println("output start");
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 1;
		final ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end,
		                                                         new IntCoordinates(512, 512));
		final TPCellIntegral<ContinuousTPShapeFunction> gg =
			new TPCellIntegralViaReferenceCell<>(1, TPCellIntegralViaReferenceCell.GRAD_GRAD);
		final ArrayList<CellIntegral<TPCell, ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		final TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1),
			                              TPRightHandSideIntegral.VALUE);
		final ArrayList<RightHandSideIntegral<TPCell, ContinuousTPShapeFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		
		final int n = grid.getShapeFunctions()
		                  .size();
		final SparseMatrix s = new SparseMatrix(n, n);
		final DenseVector d = new DenseVector(n);
		grid.writeCellIntegralsToMatrix(cellIntegrals, s);
		System.out.println("cell ints");
		grid.writeCellIntegralsToRhs(rightHandSideIntegrals, d);
		System.out.println("cell ints");
		grid.writeBoundaryValuesTo(ScalarFunction.constantFunction(0), s, d);
		System.out.println("bdr");
		
		final CartesianUpFrontSchwarz<ContinuousTPShapeFunction> schwarz =
			new CartesianUpFrontSchwarz<>(s,
			                              grid,
			                              new IntCoordinates(4, 4),
			                              3,
			                              new AdditiveSubspaceCorrection<>(1, grid),
			                              new CGSolver(1e-10));
		Vector iterate = d.mul(0);
//		for (int i = 0; i < 100; i++)
//		{
//			iterate = schwarz.getSubspaceCorrection()
//			                 .apply(schwarz, iterate, d);
//			System.out.println(s.mvMul(iterate)
//			                    .sub(d)
//			                    .euclidianNorm());
//		}
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = true;
		iterate = it.solvePCG(s, schwarz, d, 1e-7);
		
		System.out.println(it.iterations);
		
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> solut =
			new ScalarFESpaceFunction<>(
				grid.getShapeFunctionMap(), iterate);
		PlotWindow.addPlot(new ScalarPlot2D(solut, grid.generatePlotPoints(30), 30));
	}
}
