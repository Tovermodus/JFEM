import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class TaylorHoodStokes2D
{
	public static void main(String[] args)
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		builder.build();
		CoordinateVector start = CoordinateVector.fromValues(0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 1;
		TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
				TPVectorCellIntegral.GRAD_GRAD);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(gradGrad);
		MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction,QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(StokesReferenceSolution.rightHandSide(), TPVectorRightHandSideIntegral.VALUE));
		
		
		PlotWindow p = new PlotWindow();
		TaylorHoodSpace grid = new TaylorHoodSpace(start, end,
			Ints.asList(3,3));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(List.of(vv, divValue), List.of(rightHandSideIntegral));
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		grid.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
		SparseMatrix mat = grid.getSystemMatrix();
		DenseVector rhs = grid.getRhs();
		IterativeSolver i = new IterativeSolver();
		i.showProgress = false;
		Vector solution1 = i.solveCG(grid.getSystemMatrix(), grid.getRhs(), 1e-10);
		
		MixedFESpaceFunction<QkQkFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		p.addPlot(new MixedPlot2D(solut, grid.generatePlotPoints(20),20));
		TaylorHoodSpace grid2 = new TaylorHoodSpace(start, end,
			Ints.asList(3,3));
		grid2.assembleCells();
		grid2.assembleFunctions(polynomialDegree);
		grid2.initializeSystemMatrix();
		grid2.initializeRhs();
		grid2.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
		System.out.println("Cell Integrals");
		grid2.evaluateCellIntegrals(List.of(vv, divValue), List.of(rightHandSideIntegral));
		System.out.println("Face Integrals");
		grid2.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		
		System.out.println(mat.sub(grid2.getSystemMatrix()).absMaxElement());
		System.out.println(rhs.sub(grid2.getRhs()).absMaxElement());
		
		p.addPlot(new MatrixPlot(mat));
		p.addPlot(new MatrixPlot(grid2.getSystemMatrix()));
		
		
		System.out.println("solve system: " + grid2.getSystemMatrix().getRows() + "×" + grid2.getSystemMatrix().getCols());
		i = new IterativeSolver();
		i.showProgress = false;
		solution1 = i.solveCG(grid2.getSystemMatrix(), grid2.getRhs(), 1e-10);
		System.out.println("solved");
		System.out.println(grid2.getSystemMatrix().sub(grid2.getSystemMatrix().transpose()).absMaxElement());
		solut =
			new MixedFESpaceFunction<>(
				grid2.getShapeFunctions(), solution1);
		p.addPlot(new MixedPlot2D(solut, grid2.generatePlotPoints(20),20));
	}
}
