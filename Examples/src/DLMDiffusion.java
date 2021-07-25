import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.TaylorHoodSpace;
import tensorproduct.*;

import java.util.List;
import java.util.Map;

public class DLMDiffusion
{
	public static void main(String[] args)
	{
		
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		
		CoordinateVector start = CoordinateVector.fromValues(0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		CoordinateVector startImmersed = CoordinateVector.fromValues(0.25, 0.25);
		CoordinateVector endImmersed = CoordinateVector.fromValues(0.75, 0.75);
		
		int polynomialDegree = 1;
		
		ContinuousTPFESpace largeGrid = new ContinuousTPFESpace(start, end,
			Ints.asList(
				10,10), polynomialDegree);
		largeGrid.assembleCells();
		largeGrid.assembleFunctions(polynomialDegree);
		ContinuousTPFESpace immersedGrid = new ContinuousTPFESpace(startImmersed, endImmersed,
			Ints.asList(8,8), polynomialDegree);
		immersedGrid.assembleCells();
		immersedGrid.assembleFunctions(polynomialDegree);
		
		ScalarFunction rho = ScalarFunction.constantFunction(1);
		ScalarFunction rho2minrho = ScalarFunction.constantFunction(1);
		ScalarFunction f = ScalarFunction.constantFunction(2);
		ScalarFunction f2minf = ScalarFunction.constantFunction(2);
		
		TPCellIntegral<ContinuousTPShapeFunction> rhogradgrad = new TPCellIntegral<>(rho,
			TPCellIntegral.GRAD_GRAD, true);
		TPFaceIntegral<ContinuousTPShapeFunction> dirichlet =
			new TPFaceIntegral<>(ScalarFunction.constantFunction(1000),
				TPFaceIntegral.VALUE_JUMP_VALUE_JUMP, true);
		TPCellIntegral<ContinuousTPShapeFunction> rho2gradgrad = new TPCellIntegral<>(rho2minrho,
			TPCellIntegral.GRAD_GRAD, true);
		
		
		TPCellIntegral<ContinuousTPShapeFunction> lagv2 = new TPCellIntegral<>(TPCellIntegral.H1);
		
		TPRightHandSideIntegral<ContinuousTPShapeFunction> fv = new TPRightHandSideIntegral<>(f,
			TPRightHandSideIntegral.VALUE,true);
		TPRightHandSideIntegral<ContinuousTPShapeFunction> f2minfv = new TPRightHandSideIntegral<>(f2minf,
			TPRightHandSideIntegral.VALUE,true);
		
		int n = largeGrid.getShapeFunctions().size();
		int m = immersedGrid.getShapeFunctions().size();
		
		SparseMatrix A11 = new SparseMatrix(n, n);
		SparseMatrix A22 = new SparseMatrix(m,m);
		SparseMatrix A23 =  new SparseMatrix(m,m);
		SparseMatrix A13 =  new SparseMatrix(n,m);
		
		largeGrid.writeCellIntegralsToMatrix(List.of(rhogradgrad), A11);
		largeGrid.writeFaceIntegralsToMatrix(List.of(dirichlet), A11);
		
		immersedGrid.writeCellIntegralsToMatrix(List.of(rho2gradgrad), A22);
		
		immersedGrid.writeCellIntegralsToMatrix(List.of(lagv2), A23);
		
		A23.mulInPlace(-1);
		
		for (Map.Entry<Integer, ContinuousTPShapeFunction> sf: largeGrid.getShapeFunctions().entrySet())
		{
			TPRightHandSideIntegral<ContinuousTPShapeFunction> shapeFunctionOnImmersedGrid =
				new TPRightHandSideIntegral<>(sf.getValue(), TPRightHandSideIntegral.H1, false);
			DenseVector integrals = new DenseVector(m);
			immersedGrid.writeCellIntegralsToRhs(List.of(shapeFunctionOnImmersedGrid), integrals);
			A13.addSmallMatrixAt(integrals.asMatrix().transpose(), sf.getKey(),0);
		}
		
		
		SparseMatrix A =
			new SparseMatrix(n+2*m, n+2*m);
		
		SparseMatrix T =
			new SparseMatrix(n+2*m, n+2*m);
		
		A.addSmallMatrixAt(A11,0,0);
		A.addSmallMatrixAt(A22, n,n);
		A.addSmallMatrixAt(A23, n,n+m);
		A.addSmallMatrixAt(A23.transpose(), n+m,n);
		A.addSmallMatrixAt(A13, 0,n+m);
		A.addSmallMatrixAt(A13.transpose(), n+m,0);
		
		T.addSmallMatrixAt(A11.inverse(), 0, 0);
		T.addSmallMatrixAt(SparseMatrix.identity(m), n, n);
		T.addSmallMatrixAt(SparseMatrix.identity(m), n+m, n+m);
		
		DenseVector b1 = new DenseVector(n);
		DenseVector b2 = new DenseVector(m);
		
		largeGrid.writeCellIntegralsToRhs(List.of(fv), b1);
		immersedGrid.writeCellIntegralsToRhs(List.of(f2minfv), b2);
		
		DenseVector b = new DenseVector(n+2*m);
		b.addSmallVectorAt(b1, 0);
		b.addSmallVectorAt(b2, n);
		IterativeSolver i = new IterativeSolver();
		Vector solut = i.solvePGMRES(A,T,b,1e-9);//A.solve(b);
		Vector largeSolut = solut.slice(new IntCoordinates(0), new IntCoordinates(n));
		ScalarFESpaceFunction<ContinuousTPShapeFunction> solutFun =
			new ScalarFESpaceFunction<>(
				largeGrid.getShapeFunctions(), largeSolut);
		PlotWindow p = new PlotWindow();
		p.addPlot(new MatrixPlot(A));
		p.addPlot(new MatrixPlot(T));
		p.addPlot(new ScalarPlot2D(solutFun, largeGrid.generatePlotPoints(70), 70));
	}
}
