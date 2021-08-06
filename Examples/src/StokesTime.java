import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;

import java.util.List;
import java.util.Map;

public class StokesTime
{
	public static VectorFunction createBoundaryFunction(double t)
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
			public CoordinateVector value(CoordinateVector pos)
			{
				//if(Math.abs(pos.x()) <= 1e-10 || Math.abs(pos.x()-1) <= 1e-10)
					return CoordinateVector.fromValues(0.00,0);
				///if(Math.abs(pos.x()-1) <= 1e-1)
				//	return CoordinateVector.fromValues(1,1).mul(Math.sin(t)*Math.sin(t)*10*(0.25-
				//	(0.5-pos.y())*(0.5-pos.y())));
				//return new CoordinateVector(2);
			}
		};
	}
	public static MixedFunction generateCurrentFunction(Vector iterate, TaylorHoodSpace grid)
	{
		return new MixedFESpaceFunction<>(grid.getShapeFunctions(), iterate);
	}
	public static void main(String[] args)
	{
		
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		CoordinateVector start = CoordinateVector.fromValues(0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 1;
		TaylorHoodSpace grid = new TaylorHoodSpace(start, end,
			Ints.asList(11,11));
		double reynolds = 100;
		double dt = 0.3;
		int timesteps = 10;
		int nPoints = 43;
		
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction> gradGrad =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
				ScalarFunction.constantFunction(1),
				TPVectorCellIntegral.GRAD_GRAD));
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction> valueValue =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
				ScalarFunction.constantFunction(1),
				TPVectorCellIntegral.VALUE_VALUE));
		MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction,QkQkFunction> source =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>((new ScalarFunction()
				{
					@Override
					public int getDomainDimension()
					{
						return 2;
					}
					
					@Override
					public Double value(CoordinateVector pos)
					{
						return Math.exp(-10*pos.sub(CoordinateVector.fromValues(0.5,0.5)).euclidianNorm());
					}
				}).makeIsotropicVectorFunction(),
					TPVectorRightHandSideIntegral.VALUE));
		MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction,QkQkFunction> initial =
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
					public CoordinateVector value(CoordinateVector pos)
					{
						return CoordinateVector.fromValues(0,0);
					}
				},
					TPVectorRightHandSideIntegral.VALUE));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		List<CoordinateVector> points = grid.generatePlotPoints(nPoints);
		int n = grid.getShapeFunctions().size();
		
		ScalarFunction indicatorFunction = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				//if(Math.abs(pos.x()) <= 1e-10)
				//	return 1.0;
				//if(Math.abs(pos.x()-1) <= 1e-10)
				//	return 1.0;
				if(Math.abs(pos.y()) <= 1e-10)
					return 1.0;
				if(Math.abs(1-pos.y()) <= 1e-10)
					return 1.0;
				return 0.0;
			}
		};
		
		
		SparseMatrix M = new SparseMatrix(n,n);
		grid.writeCellIntegralsToMatrix(List.of(valueValue), M);
		
		SparseMatrix A = new SparseMatrix(n,n);
		grid.writeCellIntegralsToMatrix(List.of(gradGrad), A);
		A.mulInPlace(reynolds*dt);
		
		SparseMatrix D = new SparseMatrix(n,n);
		grid.writeCellIntegralsToMatrix(List.of(divValue), D);
		D.mulInPlace(100);
		
		SparseMatrix C = new SparseMatrix(n,n);
		
		DenseVector src = new DenseVector(n);
		grid.writeCellIntegralsToRhs(List.of(source), src);
		src.mulInPlace(dt);
		
		DenseVector iterate = new DenseVector(n);
		grid.writeCellIntegralsToRhs(List.of(initial), iterate);
		
		SparseMatrix MAD = M.add(A).add(D);
		
		VectorFunction bdrFunction = createBoundaryFunction(0);
		
		grid.writeBoundaryValuesTo(new MixedFunction(bdrFunction),
			(f) -> TPFaceIntegral.integrateNonTensorProduct(indicatorFunction::value, f,
				QuadratureRule1D.Gauss5) > 0, MAD, iterate);
		System.out.println(iterate);
		
		Map<CoordinateVector, Double> pvals =
			(generateCurrentFunction(iterate, grid).pressureValuesInPointsAtTime(points,0));
		Map<CoordinateVector, Double> divvals =
			(generateCurrentFunction(iterate, grid).getVelocityFunction().getDivergenceFunction().valuesInPointsAtTime(points,0));
		
		Map<CoordinateVector, CoordinateVector>vvals = (new MixedFESpaceFunction<>(
			grid.getShapeFunctions(), iterate)
			.velocityValuesInPointsAtTime(points,0));
		PlotWindow p = new PlotWindow();
		p.addPlot(new MixedPlot2D(generateCurrentFunction(iterate, grid), points, nPoints));
		for (int i = 1; i < timesteps; i++)
		{
			IterativeSolver its = new IterativeSolver();
			its.showProgress = false;
			
			
//			MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction> convection =
//				MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
//					generateCurrentFunction(iterate, grid).getVelocityFunction(),
//					TPVectorCellIntegral.GRAD_VALUE));
//			grid.writeCellIntegralsToMatrix(List.of(convection), C);
//			C.mulInPlace(dt);
			SparseMatrix MADC = MAD;//.add(C);
			
			DenseVector rhs = src.add(M.mvMul(iterate));
			bdrFunction = createBoundaryFunction(i*dt);
			grid.writeBoundaryValuesTo(new MixedFunction(bdrFunction),
				(f) -> {
					return TPFaceIntegral.integrateNonTensorProduct(indicatorFunction::value, f,
						QuadratureRule1D.Gauss5) > 0;
				}, MADC, iterate);
			iterate = new DenseVector(its.solveBiCGStab(MADC, rhs,iterate, 1e-6));
			System.out.println("x"+iterate);
			System.out.println(i);
			pvals.putAll(generateCurrentFunction(iterate, grid)
				.pressureValuesInPointsAtTime(points, dt*i));
			vvals.putAll(generateCurrentFunction(iterate, grid)
				.velocityValuesInPointsAtTime(points, dt*i));
			divvals.putAll(generateCurrentFunction(iterate, grid)
				.getVelocityFunction()
				.getDivergenceFunction().valuesInPointsAtTime(points,dt*i));
		}
		p.addPlot(new MixedPlot2DTime(pvals, vvals, nPoints));
		p.addPlot(new ScalarPlot2DTime(divvals, nPoints, ""));
		
	}
}
