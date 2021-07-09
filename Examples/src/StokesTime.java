import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.*;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;

public class StokesTime
{
	public static VectorFunction createBoundaryFunction(double t)
	{
		
		return new VectorFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				if(Math.abs(pos.x()) <= 1e-10)
					return CoordinateVector.fromValues(1,-1).mul(Math.sin(t)*Math.sin(t)*10*(0.25-(0.5-pos.y())*(0.5-pos.y())));
				///if(Math.abs(pos.x()-1) <= 1e-1)
				//	return CoordinateVector.fromValues(1,1).mul(Math.sin(t)*Math.sin(t)*10*(0.25-
				//	(0.5-pos.y())*(0.5-pos.y())));
				return new CoordinateVector(2);
			}
		};
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
			Ints.asList(6,6), polynomialDegree);
		
		MixedCellIntegral<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction, ContinuousTPVectorFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction, ContinuousTPVectorFunction> gradGrad =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<ContinuousTPVectorFunction>(
				ScalarFunction.constantFunction(1),
				TPVectorCellIntegral.GRAD_GRAD));
		MixedCellIntegral<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction, ContinuousTPVectorFunction> valueValue =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<ContinuousTPVectorFunction>(
				ScalarFunction.constantFunction(1),
				TPVectorCellIntegral.VALUE_VALUE));
		MixedRightHandSideIntegral<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction> source =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<ContinuousTPVectorFunction>(ScalarFunction.constantFunction(0).makeIsotropicVectorFunction(),
					TPVectorRightHandSideIntegral.VALUE));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		int n = grid.getShapeFunctions().size();
		double reynolds = 0.0000001;
		double dt = 0.1;
		int timesteps = 150;
		
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
				if(Math.abs(pos.x()) <= 1e-10)
					return 1.0;
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
		D.mulInPlace(dt);
		DenseVector src = new DenseVector(n);
		grid.writeCellIntegralsToRhs(List.of(source), src);
		src.mulInPlace(dt);
		SparseMatrix MAD = M.add(A).add(D);
		DenseVector iterate = new DenseVector(n);
		VectorFunction bdrFunction = createBoundaryFunction(0);
		grid.setVelocityBoundaryValues(bdrFunction, indicatorFunction, MAD);
		grid.setVelocityBoundaryValues(bdrFunction, indicatorFunction, iterate);
		
		List<CoordinateVector> points = grid.generatePlotPoints(30);
		Map<CoordinateVector, Double>pvals = (new MixedFESpaceFunction<>(
			grid.getShapeFunctions(), iterate)
			.pressureValuesInPointsAtTime(points,0));
		Map<CoordinateVector, CoordinateVector>vvals = (new MixedFESpaceFunction<>(
			grid.getShapeFunctions(), iterate)
			.velocityValuesInPointsAtTime(points,0));
		for (int i = 0; i < timesteps; i++)
		{
			IterativeSolver<Matrix> its = new IterativeSolver<>();
			its.showProgress = false;
			DenseVector rhs = src.add(M.mvMul(iterate));
			iterate = new DenseVector(its.solveBiCGStab(MAD, rhs,iterate, 1e-6));
			//System.out.println(rhs);
			grid.setVelocityBoundaryValues(bdrFunction, indicatorFunction, iterate);
			System.out.println("x"+iterate);
			System.out.println(i);
			pvals.putAll(new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), iterate)
				.pressureValuesInPointsAtTime(points, dt*i));
			vvals.putAll(new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), iterate)
				.velocityValuesInPointsAtTime(points, dt*i));
			bdrFunction = createBoundaryFunction(i*dt);
		}
		PlotWindow p = new PlotWindow();
		p.addPlot(new MixedPlot2DTime(pvals, vvals, 30));
		
	}
}
