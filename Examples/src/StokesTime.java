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
	public static void main(String[] args)
	{
		
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		CoordinateVector start = CoordinateVector.fromValues(0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 2;
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
				new TPVectorRightHandSideIntegral<ContinuousTPVectorFunction>(ScalarFunction.constantFunction(1.1).makeIsotropicVectorFunction(),
					TPVectorRightHandSideIntegral.VALUE));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		int n = grid.getShapeFunctions().size();
		double reynolds = 1;
		double dt = 0.01;
		int timesteps = 5;
		
		
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
		
		VectorFunction bdrFunction = new VectorFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				if(Math.abs(pos.x()) <= 1e-10 && pos.y() > 0.25 && pos.y() < 0.75)
					return CoordinateVector.fromValues(1,-0.5);
				if(Math.abs(pos.x()-1) <= 1e-10 && pos.y() > 0.25 && pos.y() < 0.75)
					return CoordinateVector.fromValues(1,0.9);
				return new CoordinateVector(2);
			}
		};
		grid.setVelocityBoundaryValues(bdrFunction, MAD);
		grid.setVelocityBoundaryValues(bdrFunction, iterate);
		
		List<CoordinateVector> points = grid.generatePlotPoints(30);
		Map<CoordinateVector, Double>pvals = (new MixedFESpaceFunction<>(
			grid.getShapeFunctions(), iterate)
			.pressureValuesInPointsAtTime(points,0));
		Map<CoordinateVector, Double>vvals1 = (new MixedFESpaceFunction<>(
			grid.getShapeFunctions(), iterate)
			.velocityComponentsInPointsAtTime(points,0,0));
		Map<CoordinateVector, Double>vvals2 = (new MixedFESpaceFunction<>(
			grid.getShapeFunctions(), iterate)
			.velocityComponentsInPointsAtTime(points,1,0));
		for (int i = 0; i < timesteps; i++)
		{
			IterativeSolver<Matrix> its = new IterativeSolver<>();
			its.showProgress = true;
			DenseVector rhs = src.add(M.mvMul(iterate));
			iterate = new DenseVector(its.solveGMRES(MAD, rhs, 1e-7));
			System.out.println(rhs);
			grid.setVelocityBoundaryValues(bdrFunction, iterate);
			System.out.println("x"+iterate);
			System.out.println(i);
			pvals.putAll(new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), iterate)
				.pressureValuesInPointsAtTime(points, dt*i));
			vvals1.putAll(new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), iterate)
				.velocityComponentsInPointsAtTime(points, 0,dt*i));
			vvals2.putAll(new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), iterate)
				.velocityComponentsInPointsAtTime(points, 1,dt*i));
		}
		new PlotFrame(List.of(pvals,vvals1,vvals2),start.addTime(0),end.addTime(timesteps*dt), timesteps);
		
	}
}
