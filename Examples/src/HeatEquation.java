import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.*;
import linalg.Vector;
import tensorproduct.*;

import java.util.*;

public class HeatEquation
{
	public static void main(String[] args)
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		CoordinateVector start = CoordinateVector.fromValues(-1,-1);
		CoordinateVector end = CoordinateVector.fromValues(1,1);
		int polynomialDegree = 1;
		TPFESpace grid = new TPFESpace(start,end,
			Ints.asList(12,12),polynomialDegree);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		TPCellIntegral<TPShapeFunction> gg =
			new TPCellIntegralViaReferenceCell<TPShapeFunction>(1,
				TPCellIntegral.GRAD_GRAD,
				false);
		TPCellIntegral<TPShapeFunction> vv =
			new TPCellIntegralViaReferenceCell<TPShapeFunction>(1,
				TPCellIntegral.VALUE_VALUE,
				false);
		double penalty = 200000;
		TPFaceIntegral<TPShapeFunction> jj = new TPFaceIntegral<>(ScalarFunction.constantFunction(penalty),
			TPFaceIntegral.VALUE_JUMP_VALUE_JUMP, false);
		ScalarFunction sourceFun = new ScalarFunction()
		{
			final CoordinateVector c =  CoordinateVector.fromValues(0.5,0.5);
			final CoordinateVector c2 = CoordinateVector.fromValues(-0.5,-0.5);
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return 1./(1+1e3*pos.sub(c).euclidianNorm()*pos.sub(c2).euclidianNorm());
			}
		};
		TPRightHandSideIntegral<TPShapeFunction> src =new TPRightHandSideIntegral<>(sourceFun
			, TPRightHandSideIntegral.VALUE, false);
		double dt = 0.001;
		int n = grid.getShapeFunctions().size();
		Vector iterate = new DenseVector(n);
		ScalarFESpaceFunction<TPShapeFunction> u_t;
		SparseMatrix M = new SparseMatrix(n, n);
		grid.writeCellIntegralsToMatrix(List.of(vv), M);
		grid.writeFaceIntegralsToMatrix(List.of(jj), M);
		SparseMatrix A = new SparseMatrix(n,n);
		grid.writeCellIntegralsToMatrix(List.of(gg), A);
		grid.writeFaceIntegralsToMatrix(List.of(jj), A);
		A.mulInPlace(dt);
		DenseVector source = new DenseVector(n);
		grid.writeCellIntegralsToRhs(List.of(src), source);
		source.mulInPlace(dt);
		System.out.println("source"+source);
		//System.out.println("M"+M);
		//System.out.println("A"+A);
		int timesteps = 40;
		Map<CoordinateVector, Double> vals;
		List<CoordinateVector> points = grid.generatePlotPoints(timesteps);
		vals = (new ScalarFESpaceFunction<>(
			grid.getShapeFunctions(), iterate)
			.valuesInPointsAtTime(points,0));
		System.out.println("x"+iterate);
		for (int i = 0; i < timesteps; i++)
		{
			System.out.println("Ax"+A.mvMul(iterate));
			System.out.println("Mx"+M.mvMul(iterate));
			DenseVector rhs = source.add(M.mvMul(iterate)).sub(A.mvMul(iterate));
			IterativeSolver<SparseMatrix> its = new IterativeSolver<>();
			iterate = its.solveCG(M,rhs,1e-10);
			System.out.println("x"+iterate);
			vals.putAll(new ScalarFESpaceFunction<>(
				grid.getShapeFunctions(), iterate)
				.valuesInPointsAtTime(points, dt*i));
		}
		new PlotFrame(List.of(vals),start.addTime(0),end.addTime(timesteps*dt));
		
		
	}
}
