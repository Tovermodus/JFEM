import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import tensorproduct.*;

import java.util.List;
import java.util.Map;

public class DLMDiffusion
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder
			= new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final CoordinateVector startImmersed = CoordinateVector.fromValues(0.25, 0.25);
		final CoordinateVector endImmersed = CoordinateVector.fromValues(0.75, 0.75);
		
		final int polynomialDegree = 2;
		
		final ContinuousTPFESpace largeGrid = new ContinuousTPFESpace(start, end, Ints.asList(15, 15));
		largeGrid.assembleCells();
		largeGrid.assembleFunctions(polynomialDegree);
		final ContinuousTPFESpace immersedGrid = new ContinuousTPFESpace(startImmersed, endImmersed,
		                                                                 Ints.asList(8, 8));
		immersedGrid.assembleCells();
		immersedGrid.assembleFunctions(polynomialDegree);
		
		final ScalarFunction rho = ScalarFunction.constantFunction(1);
		final ScalarFunction rho2minrho = ScalarFunction.constantFunction(100);
		final ScalarFunction f = ScalarFunction.constantFunction(2);
		final ScalarFunction f2minf = ScalarFunction.constantFunction(-4);
		
		final TPCellIntegral<ContinuousTPShapeFunction> rhogradgrad = new TPCellIntegral<>(rho,
		                                                                                   TPCellIntegral.GRAD_GRAD);
		final TPFaceIntegral<ContinuousTPShapeFunction> dirichlet = new TPFaceIntegral<>(
			ScalarFunction.constantFunction(1000), TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
		final TPCellIntegral<ContinuousTPShapeFunction> rho2gradgrad = new TPCellIntegral<>(rho2minrho,
		                                                                                    TPCellIntegral.GRAD_GRAD);
		
		final TPCellIntegral<ContinuousTPShapeFunction> lagv2 = new TPCellIntegral<>(TPCellIntegral.H1);
		
		final TPRightHandSideIntegral<ContinuousTPShapeFunction> fv = new TPRightHandSideIntegral<>(f,
		                                                                                            TPRightHandSideIntegral.VALUE);
		final TPRightHandSideIntegral<ContinuousTPShapeFunction> f2minfv = new TPRightHandSideIntegral<>(f2minf,
		                                                                                                 TPRightHandSideIntegral.VALUE);
		
		final int n = largeGrid.getShapeFunctions().size();
		final int m = immersedGrid.getShapeFunctions().size();
		
		final SparseMatrix A11 = new SparseMatrix(n, n);
		final SparseMatrix A22 = new SparseMatrix(m, m);
		final SparseMatrix A23 = new SparseMatrix(m, m);
		final SparseMatrix A13 = new SparseMatrix(n, m);
		
		largeGrid.writeCellIntegralsToMatrix(List.of(rhogradgrad), A11);
		largeGrid.writeFaceIntegralsToMatrix(List.of(dirichlet), A11);
		
		immersedGrid.writeCellIntegralsToMatrix(List.of(rho2gradgrad), A22);
		
		immersedGrid.writeCellIntegralsToMatrix(List.of(lagv2), A23);
		
		A23.mulInPlace(-1);
		
		for (final Map.Entry<Integer, ContinuousTPShapeFunction> sf : largeGrid.getShapeFunctions().entrySet())
		{
			final TPRightHandSideIntegral<ContinuousTPShapeFunction> shapeFunctionOnImmersedGrid
				= new TPRightHandSideIntegral<>(sf.getValue(), TPRightHandSideIntegral.H1);
			final DenseVector integrals = new DenseVector(m);
			immersedGrid.writeCellIntegralsToRhs(List.of(shapeFunctionOnImmersedGrid), integrals);
			A13.addSmallMatrixAt(integrals.asMatrix().transpose(), sf.getKey(), 0);
		}
		
		final SparseMatrix A = new SparseMatrix(n + 2 * m, n + 2 * m);
		
		final SparseMatrix T = new SparseMatrix(n + 2 * m, n + 2 * m);
		
		A.addSmallMatrixAt(A11, 0, 0);
		A.addSmallMatrixAt(A22, n, n);
		A.addSmallMatrixAt(A23, n, n + m);
		A.addSmallMatrixAt(A23.transpose(), n + m, n);
		A.addSmallMatrixAt(A13, 0, n + m);
		A.addSmallMatrixAt(A13.transpose(), n + m, 0);
		
		final DirectlySolvable A11inv = A11.inverse();
		T.addSmallMatrixAt(A11inv, 0, 0);
		T.addSmallMatrixAt((A22.add(SparseMatrix.identity(m).mul(0.01))).inverse(), n, n);
		T.addSmallMatrixAt(SparseMatrix.identity(m), n + m, n + m);
		
		final DenseVector b1 = new DenseVector(n);
		final DenseVector b2 = new DenseVector(m);
		
		largeGrid.writeCellIntegralsToRhs(List.of(fv), b1);
		immersedGrid.writeCellIntegralsToRhs(List.of(f2minfv), b2);
		
		final DenseVector b = new DenseVector(n + 2 * m);
		b.addSmallVectorAt(b1, 0);
		b.addSmallVectorAt(b2, n);
		final IterativeSolver i = new IterativeSolver();
		final Vector solut = i.solvePGMRES(A, T, b, 1e-9);//A.solve(b);
		final Vector largeSolut = solut.slice(new IntCoordinates(0), new IntCoordinates(n));
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> solutFun = new ScalarFESpaceFunction<>(
			largeGrid.getShapeFunctions(), largeSolut);
		final PlotWindow p = new PlotWindow();
		p.addPlot(new ScalarPlot2D(solutFun, largeGrid.generatePlotPoints(70), 70));
	}
}
