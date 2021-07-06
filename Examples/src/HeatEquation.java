import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.SparseMatrix;
import linalg.Vector;
import tensorproduct.*;

import java.util.ArrayList;

public class HeatEquation
{
	public static void main(String[] args)
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		CoordinateVector start = CoordinateVector.fromValues(-1,-1);
		CoordinateVector end = CoordinateVector.fromValues(1,1);
		int polynomialDegree = 2;
		TPFESpace grid = new TPFESpace(start,end,
			Ints.asList(5,5),polynomialDegree);
		TPCellIntegral<TPShapeFunction> gg =
			new TPCellIntegralViaReferenceCell<TPShapeFunction>(1,
			TPCellIntegral.GRAD_GRAD,
			false);
		double penalty = 200000;
		TPFaceIntegral<TPShapeFunction> jj = new TPFaceIntegral<>(ScalarFunction.constantFunction(penalty),
			TPFaceIntegral.VALUE_JUMP_VALUE_JUMP, false);
		ArrayList<CellIntegral<TPCell,TPShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		ArrayList<FaceIntegral<TPFace,TPShapeFunction>> faceIntegrals = new ArrayList<>();
		faceIntegrals.add(jj);
		TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(0),TPRightHandSideIntegral.VALUE,
				true);
		ArrayList<RightHandSideIntegral<TPCell,TPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		TPBoundaryFaceIntegral<TPShapeFunction> bound = new TPBoundaryFaceIntegral<>(new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				if (Math.abs(pos.x()) == 1||Math.abs(pos.y()) == 1)
					return +penalty*2*(1+pos.y())/((3+pos.x())*(3+pos.x())+(1+pos.y())*(1+pos.y()));
				return (double) 0;
			}
		},TPBoundaryFaceIntegral.VALUE,false);
		
		ArrayList<BoundaryRightHandSideIntegral<TPFace,TPShapeFunction>> boundaryFaceIntegrals = new ArrayList<>();
		boundaryFaceIntegrals.add(bound);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals);
		System.out.println("solve system: "+grid.getSystemMatrix().getRows()+"×"+grid.getSystemMatrix().getCols());
		//grid.A.makeParallelReady(12);
		if(grid.getRhs().getLength() < 50)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		IterativeSolver<SparseMatrix> i = new IterativeSolver<>();
		System.out.println("start stopwatch");
		Stopwatch s = Stopwatch.createStarted();
		Vector solution1 = i.solveCG(grid.getSystemMatrix(),grid.getRhs(),1e-3);
		
	}
}
