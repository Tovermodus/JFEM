import basic.PerformanceArguments;
import basic.PlotWindow;
import basic.ScalarFunction;
import basic.VectorFunction;
import com.google.common.primitives.Ints;
import distorted.CircleOverlay;
import distorted.DistortedVectorCellIntegral;
import distorted.DistortedVectorRightHandSideIntegral;
import distorted.DistortedVectorSpace;
import linalg.*;
import mixed.*;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.geometry.TPCell;

import java.util.List;
import java.util.Map;

public class DLMStokesTime
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder
			= new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		
		final CoordinateVector start = CoordinateVector.fromValues(0, -2);
		final CoordinateVector end = CoordinateVector.fromValues(5, 7);
		final CoordinateVector immercedCenter = CoordinateVector.fromValues(1.3, 2.4);
		final double immersedRadius = 0.8;
		final double reynoldsNumber = 1000;
		final int polynomialDegree = 1;
		final int timesteps = 100;
		final double dt = 0.1;
		
		final TaylorHoodSpace largeGrid = new TaylorHoodSpace(start, end, Ints.asList(12, 12));
		
		final int nPoints = 40;
		final List<CoordinateVector> points = largeGrid.generatePlotPoints(nPoints);
		largeGrid.assembleCells();
		largeGrid.assembleFunctions(polynomialDegree);
		final DistortedVectorSpace immersedGrid = new DistortedVectorSpace(immercedCenter, immersedRadius, 0);
		immersedGrid.assembleCells();
		immersedGrid.assembleFunctions(polynomialDegree);
		
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			reynolds = MixedCellIntegral.fromVelocityIntegral(
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(reynoldsNumber),
			                           TPVectorCellIntegral.SYM_GRAD));
		final MixedTPCellIntegral<ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> divValue
			= new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1), MixedTPCellIntegral.DIV_VALUE);
		
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
		final MixedFunction boundaryValues = new MixedFunction(new VectorFunction()
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
				return CoordinateVector.fromValues(1000, 0);
			}
		});
		largeGrid.writeBoundaryValuesTo(boundaryValues, (f) -> f.center().x() == 0,
		                                (f, sf) -> sf.hasVelocityFunction(), A11, b1);
		final MixedFunction pboundaryValues = new MixedFunction(new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return 0.;
			}
		});
		largeGrid.writeBoundaryValuesTo(pboundaryValues, (f) -> f.center().x() == 0 && f.center().y() > 4.53,
		                                (f, sf) -> sf.hasPressureFunction(), A11, b1);
		System.out.println("A11");
		A22 = SparseMatrix.identity(m);
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
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return CoordinateDenseMatrix.fromValues(2, 2, 1, 0, 0, 1);
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
					return sf
						.getValue()
						.getVelocityShapeFunction()
						.value(elasticTransformation.value(pos));
				}
				
				@Override
				public CoordinateMatrix gradient(final CoordinateVector pos)
				{
					return sf
						.getValue()
						.getVelocityShapeFunction()
						.gradient(
							elasticTransformation.value(pos));//.mvmul(elastictransformation
					// .gradient oder so
				}
			};
			if (sf.getValue().hasVelocityFunction())
			{
				final DistortedVectorRightHandSideIntegral shapeFunctionOnImmersedGrid
					= new DistortedVectorRightHandSideIntegral(toBeMultiplierd,
					                                           DistortedVectorRightHandSideIntegral.H1);
				final DenseVector integrals = new DenseVector(m);
				immersedGrid.writeCellIntegralsToRhs(List.of(shapeFunctionOnImmersedGrid), integrals);
				A13.addSmallMatrixAt(integrals.asMatrix().transpose(), sf.getKey(), 0);
				//if (count % 10 == 0)
				System.out.println("A31: " + 100.0 * (count++) / n + "%");
			}
		}
		
		final SparseMatrix A = new SparseMatrix(n + 2 * m, n + 2 * m);
		
		final SparseMatrix T = new SparseMatrix(n + 2 * m, n + 2 * m);
		
		A.addSmallMatrixAt(A11, 0, 0);
		A.addSmallMatrixAt(A22, n, n);
		//A.addSmallMatrixAt(A23, n, n + m);
		A.addSmallMatrixAt(A23.transpose(), n + m, n);
		A.addSmallMatrixAt(A13, 0, n + m);
		A.addSmallMatrixAt(A13.transpose(), n + m, 0);
		for (final int i : largeGrid.getFixedNodeIndices())
			A.deleteRow(i);
		for (final int i : immersedGrid.getFixedNodeIndices())
			A.deleteRow(i + n);
		System.out.println("A");
		//final DenseMatrix A11inv = A11.inverse();
		//T.addSmallMatrixAt(A11inv, 0, 0);
		System.out.println("T11");
		//T.addSmallMatrixAt((A22.add(SparseMatrix.identity(m).mul(0.1))).inverse(), n, n);
		System.out.println("T22");
		//T.addSmallMatrixAt(SparseMatrix.identity(m), n + m, n + m);
		System.out.println("T");
		
		final DenseVector b = new DenseVector(n + 2 * m);
		b.addSmallVectorAt(b1, 0);
		b.addSmallVectorAt(b2, n);
		System.out.println("b");
		
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			valueValue = MixedCellIntegral.fromVelocityIntegral(
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(1),
			                           TPVectorCellIntegral.VALUE_VALUE));
		
		final SparseMatrix M11 = new SparseMatrix(n, n);
		largeGrid.writeCellIntegralsToMatrix(List.of(valueValue), M11);
		final SparseMatrix M = new SparseMatrix(n + 2 * m, n + 2 * m);
		M.addSmallMatrixAt(M11, 0, 0);
		
		A.mulInPlace(dt);
		b.mulInPlace(dt);
		A.addInPlace(M);
		
		Vector iterate = new DenseVector(n + 2 * m);
		
		final Map<CoordinateVector, Double> pvals = (new MixedFESpaceFunction<>(largeGrid.getShapeFunctions(),
		                                                                        iterate).pressureValuesInPointsAtTime(
			points, 0));
		final Map<CoordinateVector, CoordinateVector> vvals = (new MixedFESpaceFunction<>(
			largeGrid.getShapeFunctions(), iterate).velocityValuesInPointsAtTime(points, 0));
		System.out.println("inverting");
		final DirectlySolvable Ainv = A.inverse();
		System.out.println("inverted");
		for (int i = 1; i < timesteps; i++)
		{
			System.out.println(i);
			final Vector rhs = b.add(M.mvMul(iterate));
			iterate = Ainv.mvMul(rhs);
			vvals.putAll((new MixedFESpaceFunction<>(largeGrid.getShapeFunctions(),
			                                         iterate).velocityValuesInPointsAtTime(points,
			                                                                               i * dt)));
			pvals.putAll((new MixedFESpaceFunction<>(largeGrid.getShapeFunctions(),
			                                         iterate).pressureValuesInPointsAtTime(points,
			                                                                               i * dt)));
		}
		final PlotWindow p = new PlotWindow();
		final MixedPlot2DTime plot = new MixedPlot2DTime(pvals, vvals, nPoints);
		p.addPlot(plot.addOverlay(new CircleOverlay(immersedGrid, plot)));
	}
}
