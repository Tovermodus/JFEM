import basic.*;
import com.google.common.base.Stopwatch;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.IterativeSolver;
import linalg.Vector;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeSet;

public class LaplaceContinuous
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		
		final long startTime = System.nanoTime();
		
		System.out.println("output start");
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 3;
		final ContinuousTPFESpace grid = new ContinuousTPFESpace(start, end,
		                                                         new IntCoordinates(8, 8));
		final TPCellIntegral<ContinuousTPShapeFunction> gg =
			new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			                     TPCellIntegral.GRAD_GRAD);
		final double penalty = 1;
		final TPFaceIntegral<ContinuousTPShapeFunction> jj =
			new TPFaceIntegral<>(ScalarFunction.constantFunction(penalty),
			                     TPFaceIntegral.BOUNDARY_VALUE);
		final ArrayList<CellIntegral<TPCell, ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		final ArrayList<FaceIntegral<TPFace, ContinuousTPShapeFunction>> faceIntegrals = new ArrayList<>();
		faceIntegrals.add(jj);
		final TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(4),
			                              TPRightHandSideIntegral.VALUE);
		final ArrayList<RightHandSideIntegral<TPCell, ContinuousTPShapeFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		final ArrayList<BoundaryRightHandSideIntegral<TPFace, ContinuousTPShapeFunction>> boundaryFaceIntegrals
			= new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		//grid.getSystemMatrix().add(1,0,0);
		System.out.println(((1.0 * System.nanoTime() - startTime) / 1e9));
		System.out.println("solve system: " + grid.getSystemMatrix()
		                                          .getRows() + "×" + grid.getSystemMatrix()
		                                                                 .getCols());
		//grid.A.makeParallelReady(12);
		if (grid.getRhs()
		        .getLength() < 50)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		final IterativeSolver it = new IterativeSolver();
		System.out.println("start stopwatch");
		final Stopwatch s = Stopwatch.createStarted();
		final Vector solution1 = it.solveCG(grid.getSystemMatrix(), grid.getRhs(), 1e-3);
		System.out.println(s.elapsed());
		//Vector solution = ((DenseMatrix)grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		System.out.println(((1.0 * System.nanoTime() - startTime) / 1e9));
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> solut =
			new ScalarFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		final Map<CoordinateVector, Double> vals = solut.valuesInPoints(grid.generatePlotPoints(50));
		final ArrayList<Map<CoordinateVector, Double>> valList = new ArrayList<>();
		final TreeSet<ContinuousTPShapeFunction> shapeFunctionTreeSet =
			new TreeSet<>(grid.getShapeFunctions()
			                  .values());
		valList.add(solut.valuesInPoints(grid.generatePlotPoints(50)));
		for (final ContinuousTPShapeFunction sf : shapeFunctionTreeSet)
			valList.add(sf.valuesInPoints(grid.generatePlotPoints(50)));
		new PlotFrame(valList, start, end);
	}
}
