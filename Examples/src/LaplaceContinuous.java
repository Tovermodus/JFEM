import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.*;
import tensorproduct.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.ArrayList;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.TreeSet;

public class LaplaceContinuous
{
	public static void main(String[] args)
	{
		
	
//[+6.667e+02  +1.667e+02  +1.667e+02  -3.333e-03
//		+1.667e+02  +6.667e+02  -3.333e-03  +1.667e+02
//		+1.667e+02  -3.333e-03  +6.667e+02  +1.667e+02
//		-3.333e-03  +1.667e+02  +1.667e+02  +6.667e+02
// ]
		long startTime = System.nanoTime();
		
		System.out.println("output start");
		CoordinateVector start = CoordinateVector.fromValues(-1,-1);
		CoordinateVector end = CoordinateVector.fromValues(1,1);
		int polynomialDegree = 3;
		ContinuousTPFESpace grid = new ContinuousTPFESpace(start,end,
			Ints.asList(3,3),polynomialDegree);
		TPCellIntegral<ContinuousTPShapeFunction> gg =
			new TPCellIntegral<>(ScalarFunction.constantFunction(1),
			TPCellIntegral.GRAD_GRAD,
			false);
		double penalty = 1;
		TPFaceIntegral<ContinuousTPShapeFunction> jj =
			new TPFaceIntegral<>(ScalarFunction.constantFunction(penalty),
			TPFaceIntegral.BOUNDARY_VALUE, false);
		ArrayList<CellIntegral<TPCell,ContinuousTPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		ArrayList<FaceIntegral<TPFace,ContinuousTPShapeFunction>> faceIntegrals = new ArrayList<>();
		faceIntegrals.add(jj);
		TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(4),TPRightHandSideIntegral.VALUE,
				true);
		ArrayList<RightHandSideIntegral<TPCell,ContinuousTPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		ArrayList<BoundaryRightHandSideIntegral<TPFace,ContinuousTPShapeFunction>> boundaryFaceIntegrals = new ArrayList<>();
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
		IterativeSolver<SparseMatrix> it = new IterativeSolver<>();
		System.out.println("start stopwatch");
		Stopwatch s = Stopwatch.createStarted();
		Vector solution1 = it.solveCG(grid.getSystemMatrix(),grid.getRhs(),1e-3);
		System.out.println(s.elapsed());
		//Vector solution = ((DenseMatrix)grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		ScalarFESpaceFunction<ContinuousTPShapeFunction> solut =
			new ScalarFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		Map<CoordinateVector, Double> vals = solut.valuesInPoints(grid.generatePlotPoints(50));
		ArrayList<Map<CoordinateVector, Double>> valList = new ArrayList<>();
		TreeSet<ContinuousTPShapeFunction> shapeFunctionTreeSet =
			new TreeSet<>(grid.getShapeFunctions().values());
		valList.add(solut.valuesInPoints(grid.generatePlotPoints(50)));
		for(ContinuousTPShapeFunction sf:shapeFunctionTreeSet)
			valList.add(sf.valuesInPoints(grid.generatePlotPoints(50)));
		new PlotFrame(valList, start, end);
	}
}
