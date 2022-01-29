import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class LaplaceScalarRT
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
		final int polynomialDegree = 4;
		final ScalarRTFESpace grid = new ScalarRTFESpace(start, end,
		                                                 Ints.asList(7, 7));
		final TPCellIntegral<RTComponentFunction> gg =
			new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			                     TPCellIntegral.GRAD_GRAD);
		final TPFaceIntegral<RTComponentFunction> jj = new TPFaceIntegral<>(ScalarFunction.constantFunction(10),
		                                                                    TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
		final ArrayList<CellIntegral<TPCell, RTComponentFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		final ArrayList<FaceIntegral<TPFace, RTComponentFunction>> faceIntegrals = new ArrayList<>();
		faceIntegrals.add(jj);
		final TPRightHandSideIntegral<RTComponentFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(4),
			                              TPRightHandSideIntegral.VALUE);
		final ArrayList<RightHandSideIntegral<TPCell, RTComponentFunction>> rightHandSideIntegrals
			= new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		final ArrayList<BoundaryRightHandSideIntegral<TPFace, RTComponentFunction>> boundaryFaceIntegrals
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
		final ScalarFESpaceFunction<RTComponentFunction> solut =
			new ScalarFESpaceFunction<>(
				grid.getShapeFunctionMap(), solution1);
		final Map<String, Map<CoordinateVector, Double>> valList = new HashMap<>();
		final TreeSet<RTComponentFunction> shapeFunctionTreeSet =
			new TreeSet<>(grid.getShapeFunctionMap()
			                  .values());
		valList.put("solution", solut.valuesInPoints(grid.generatePlotPoints(50)));
		for (final RTComponentFunction sf : shapeFunctionTreeSet)
			valList.put("Shapefunction", sf.valuesInPoints(grid.generatePlotPoints(50)));
		new PlotFrame(valList, start, end);
	}
}
