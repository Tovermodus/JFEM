import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.SparseMatrix;
import linalg.Vector;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class LaplaceScalarRT
{
	public static void main(String[] args)
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		long startTime = System.nanoTime();
		
		System.out.println("output start");
		CoordinateVector start = CoordinateVector.fromValues(0,0);
		CoordinateVector end = CoordinateVector.fromValues(1,1);
		int polynomialDegree = 4;
		ScalarRTFESpace grid = new ScalarRTFESpace(start,end,
			Ints.asList(7,7),polynomialDegree);
		TPCellIntegral<RTComponentFunction> gg =
			new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD);
		TPFaceIntegral<RTComponentFunction> jj = new TPFaceIntegral<>(ScalarFunction.constantFunction(10),
			TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
		ArrayList<CellIntegral<TPCell,RTComponentFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		ArrayList<FaceIntegral<TPFace,RTComponentFunction>> faceIntegrals = new ArrayList<>();
		faceIntegrals.add(jj);
		TPRightHandSideIntegral<RTComponentFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(4),TPRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell,RTComponentFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		ArrayList<BoundaryRightHandSideIntegral<TPFace,RTComponentFunction>> boundaryFaceIntegrals = new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
		grid.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals);
		//grid.getSystemMatrix().add(1,0,0);
		System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
		System.out.println("solve system: "+grid.getSystemMatrix().getRows()+"×"+grid.getSystemMatrix().getCols());
		//grid.A.makeParallelReady(12);
		if(grid.getRhs().getLength() < 50)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		IterativeSolver it = new IterativeSolver();
		System.out.println("start stopwatch");
		Stopwatch s = Stopwatch.createStarted();
		Vector solution1 = it.solveCG(grid.getSystemMatrix(),grid.getRhs(),1e-3);
		System.out.println(s.elapsed());
		//Vector solution = ((DenseMatrix)grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		ScalarFESpaceFunction<RTComponentFunction> solut =
			new ScalarFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		Map<String,Map<CoordinateVector, Double>> valList = new HashMap<>();
		TreeSet<RTComponentFunction> shapeFunctionTreeSet =
			new TreeSet<>(grid.getShapeFunctions().values());
		valList.put("solution",solut.valuesInPoints(grid.generatePlotPoints(50)));
		for(RTComponentFunction sf:shapeFunctionTreeSet)
			valList.put("Shapefunction",sf.valuesInPoints(grid.generatePlotPoints(50)));
		new PlotFrame(valList, start, end);
	}
}
