import basic.PerformanceArguments;
import basic.ScalarFunction;
import basic.VectorFunction;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;

import java.util.List;
import java.util.Map;

public class StokesTime
{
	public static VectorFunction createBoundaryFunction(final double t)
	{
		
		return new VectorFunction()
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
				if (Math.abs(pos.x()) <= 1e-10 || Math.abs(pos.x() - 1) <= 1e-10)
					return CoordinateVector.fromValues(0.1, 0);
				///if(Math.abs(pos.x()-1) <= 1e-1)
				//	return CoordinateVector.fromValues(1,1).mul(Math.sin(t)*Math.sin(t)*10*(0.25-
				//	(0.5-pos.y())*(0.5-pos.y())));
				return new CoordinateVector(2);
			}
		};
	}
	
	public static MixedFunction generateCurrentFunction(final Vector iterate, final TaylorHoodSpace grid)
	{
		return new MixedFESpaceFunction<>(grid.getShapeFunctions(), iterate);
	}
	
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 1;
		final TaylorHoodSpace grid = new TaylorHoodSpace(start, end,
		                                                 Ints.asList(5, 5));
		final double reynolds = 1;
		final double dt = 0.003;
		final int timesteps = 3;
		final int nPoints = 43;
		
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			gradGrad =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
				ScalarFunction.constantFunction(1),
				TPVectorCellIntegral.GRAD_GRAD));
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			valueValue =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
				ScalarFunction.constantFunction(1),
				TPVectorCellIntegral.VALUE_VALUE));
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> source =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(
					new VectorFunction()
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
							return CoordinateVector.fromValues(0.1,
							                                   10 * Math.sin(pos.x() * 20));
						}
					},
					TPVectorRightHandSideIntegral.VALUE));
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> initial =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(new VectorFunction()
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
						return CoordinateVector.fromValues(0, 0);
					}
				},
				                                    TPVectorRightHandSideIntegral.VALUE));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		final List<CoordinateVector> points = grid.generatePlotPoints(nPoints);
		final int n = grid.getShapeFunctions()
		                  .size();
		
		final ScalarFunction indicatorFunction = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
//				if (Math.abs(pos.x()) <= 1e-10)
//					return 1.0;
				//if(Math.abs(pos.x()-1) <= 1e-10)
				//	return 1.0;
				if (Math.abs(pos.y()) <= 1e-10)
					return 1.0;
				if (Math.abs(1 - pos.y()) <= 1e-10)
					return 1.0;
				return 0.0;
			}
		};
		
		final SparseMatrix M = new SparseMatrix(n, n);
		grid.writeCellIntegralsToMatrix(List.of(valueValue), M);
		
		final SparseMatrix A = new SparseMatrix(n, n);
		grid.writeCellIntegralsToMatrix(List.of(gradGrad), A);
		A.mulInPlace(reynolds * dt);
		
		final SparseMatrix D = new SparseMatrix(n, n);
		grid.writeCellIntegralsToMatrix(List.of(divValue), D);
		D.mulInPlace(100);
		
		final SparseMatrix C = new SparseMatrix(n, n);
		
		final DenseVector src = new DenseVector(n);
		grid.writeCellIntegralsToRhs(List.of(source), src);
		src.mulInPlace(dt);
		
		DenseVector iterate = new DenseVector(n);
		grid.writeCellIntegralsToRhs(List.of(initial), iterate);
		
		final SparseMatrix MAD = M.add(A)
		                          .add(D);
		
		VectorFunction bdrFunction = createBoundaryFunction(0);
		
		grid.writeBoundaryValuesTo(new ComposedMixedFunction(bdrFunction),
		                           (f) -> TPFaceIntegral.integrateNonTensorProduct(indicatorFunction::value, f,
		                                                                           QuadratureRule1D.Gauss5) > 0,
		                           MAD, iterate);
		System.out.println(iterate);
		
		final Map<CoordinateVector, Double> pvals =
			(generateCurrentFunction(iterate, grid).pressureValuesInPointsAtTime(points, 0));
		final Map<CoordinateVector, Double> divvals =
			(generateCurrentFunction(iterate, grid)
				 .getVelocityFunction()
				 .getDivergenceFunction()
				 .valuesInPointsAtTime(points, 0));
		
		final Map<CoordinateVector, CoordinateVector> vvals = (new MixedFESpaceFunction<>(
			grid.getShapeFunctions(), iterate)
			                                                       .velocityValuesInPointsAtTime(points,
			                                                                                     0));
		//final PlotWindow p = new PlotWindow();
		//p.addPlot(new MixedPlot2D(generateCurrentFunction(iterate, grid), points, nPoints));
		for (int i = 1; i < timesteps; i++)
		{
			final IterativeSolver its = new IterativeSolver();
			its.showProgress = false;
			
			final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
				convection =
				MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
					generateCurrentFunction(iterate, grid).getVelocityFunction(),
					TPVectorCellIntegral.GRAD_VALUE));
			grid.writeCellIntegralsToMatrix(List.of(convection), C);
			C.mulInPlace(dt);
			final SparseMatrix MADC = MAD.add(C);
			
			final DenseVector rhs = src.add(M.mvMul(iterate));
			bdrFunction = createBoundaryFunction(i * dt);
			grid.writeBoundaryValuesTo(new ComposedMixedFunction(bdrFunction),
			                           (f) ->
			                           {
				                           return TPFaceIntegral.integrateNonTensorProduct(
					                           indicatorFunction::value, f,
					                           QuadratureRule1D.Gauss5) > 0;
			                           }, MADC, iterate);
			iterate = new DenseVector(its.solveBiCGStab(MADC, rhs, iterate, 1e-6));
			System.out.println("x" + iterate);
			System.out.println(i);
//			pvals.putAll(generateCurrentFunction(iterate, grid)
//				             .pressureValuesInPointsAtTime(points, dt * i));
//			vvals.putAll(generateCurrentFunction(iterate, grid)
//				             .velocityValuesInPointsAtTime(points, dt * i));
//			divvals.putAll(generateCurrentFunction(iterate, grid)
//				               .getVelocityFunction()
//				               .getDivergenceFunction().valuesInPointsAtTime(points, dt * i));
		}
//		p.addPlot(new MixedPlot2DTime(pvals, vvals, nPoints));
//		p.addPlot(new ScalarPlot2DTime(divvals, nPoints, ""));
	}
}
