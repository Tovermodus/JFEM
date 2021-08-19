import basic.PerformanceArguments;
import basic.PlotWindow;
import basic.ScalarFunction;
import basic.VectorFunction;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.*;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.geometry.TPCell;

import java.util.List;
import java.util.Map;

public class DLMStokes
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(5, 1);
		final CoordinateVector immercedCenter = CoordinateVector.fromValues(0.5, 0.5);
		final double immersedRadius = 0.1;
		final double reynoldsNumber = 250;
		final int polynomialDegree = 2;
		
		final TaylorHoodSpace largeGrid = new TaylorHoodSpace(start, end,
		                                                      Ints.asList(
			                                                      15, 15));
		largeGrid.assembleCells();
		largeGrid.assembleFunctions(polynomialDegree);
		final DistortedVectorSpace immersedGrid = new DistortedVectorSpace(immercedCenter, immersedRadius, 2);
		immersedGrid.assembleCells();
		immersedGrid.assembleFunctions(polynomialDegree);
		
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			reynolds =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
				ScalarFunction.constantFunction(reynoldsNumber), TPVectorCellIntegral.GRAD_GRAD));
		final MixedTPCellIntegral<ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue = new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
			                                     MixedTPCellIntegral.DIV_VALUE);
		
		final DistortedVectorCellIntegral lagv2 = new DistortedVectorCellIntegral(
			DistortedVectorCellIntegral.H1);
		
		final int n = largeGrid.getShapeFunctions().size();
		final int m = immersedGrid.getShapeFunctions().size();
		
		final SparseMatrix A11 = new SparseMatrix(n, n);
		SparseMatrix A22 = new SparseMatrix(m, m);
		final SparseMatrix A23 = new SparseMatrix(m, m);
		final SparseMatrix A13 = new SparseMatrix(n, m);
		
		final DenseVector b1 = new DenseVector(n);
		final DenseVector b2 = new DenseVector(m);
		largeGrid.writeCellIntegralsToMatrix(List.of(reynolds), A11);
		largeGrid.writeCellIntegralsToMatrix(List.of(divValue), A11);
		largeGrid.writeBoundaryValuesTo(new MixedFunction(new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(1, 0);
			}
		}), (f) -> f.center().x() == 0, (f, sf) -> sf.hasVelocityFunction(), A11, b1);
		System.out.println("A11");
		A22 = SparseMatrix.identity(m);
		//immersedGrid.writeCellIntegralsToMatrix(List.of(rho2gradgrad), A22);
		System.out.println("A22");
		
		immersedGrid.writeCellIntegralsToMatrix(List.of(lagv2), A23);
		
		A23.mulInPlace(-1);
		System.out.println("A23");
		int count = 0;
		final VectorFunction elasticTransformation = new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return pos;
			}
		};
		for (final Map.Entry<Integer, QkQkFunction> sf : largeGrid.getShapeFunctions().entrySet())
		{
			final VectorFunction toBeMultiplierd = new VectorFunction() //THIS NEEDS TO BE UPDATED IN
				// EVERY ITERATION
			{
				@Override
				public int getRangeDimension()
				{
					return 2;
				}
				
				@Override
				public int getDomainDimension()
				{
					return 2;
				}
				
				@Override
				public CoordinateVector value(final CoordinateVector pos)
				{
					return sf.getValue().getVelocityShapeFunction().value(
						elasticTransformation.value(pos));
				}
			};
			final DistortedVectorRightHandSideIntegral shapeFunctionOnImmersedGrid =
				new DistortedVectorRightHandSideIntegral(toBeMultiplierd,
				                                         DistortedVectorRightHandSideIntegral.H1);
			final DenseVector integrals = new DenseVector(m);
			immersedGrid.writeCellIntegralsToRhs(List.of(shapeFunctionOnImmersedGrid), integrals);
			A13.addSmallMatrixAt(integrals.asMatrix().transpose(), sf.getKey(), 0);
			//if (count % 10 == 0)
			System.out.println("A31: " + 100.0 * (count++) / n + "%");
		}
		
		final SparseMatrix A =
			new SparseMatrix(n + 2 * m, n + 2 * m);
		
		final SparseMatrix T =
			new SparseMatrix(n + 2 * m, n + 2 * m);
		
		A.addSmallMatrixAt(A11, 0, 0);
		A.addSmallMatrixAt(A22, n, n);
		A.addSmallMatrixAt(A23, n, n + m);
		A.addSmallMatrixAt(A23.transpose(), n + m, n);
		A.addSmallMatrixAt(A13, 0, n + m);
		A.addSmallMatrixAt(A13.transpose(), n + m, 0);
		System.out.println("A");
		final DenseMatrix A11inv = A11.inverse();
		T.addSmallMatrixAt(A11inv, 0, 0);
		System.out.println("T11");
		T.addSmallMatrixAt((A22.add(SparseMatrix.identity(m).mul(0.01))).inverse(), n, n);
		System.out.println("T22");
		T.addSmallMatrixAt(SparseMatrix.identity(m), n + m, n + m);
		System.out.println("T");
		
		final DenseVector b = new DenseVector(n + 2 * m);
		b.addSmallVectorAt(b1, 0);
		b.addSmallVectorAt(b2, n);
		System.out.println("b");
		final IterativeSolver i = new IterativeSolver();
		final Vector solut = i.solvePGMRES(A, T, b, 1e-9);//A.solve(b);
		final Vector largeSolut = solut.slice(new IntCoordinates(0), new IntCoordinates(n));
		final MixedFESpaceFunction<QkQkFunction> solutFun =
			new MixedFESpaceFunction<>(
				largeGrid.getShapeFunctions(), largeSolut);
		final PlotWindow p = new PlotWindow();
		p.addPlot(new MixedPlot2D(solutFun, largeGrid.generatePlotPoints(70), 70));
	}
}
