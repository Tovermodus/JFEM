import basic.VectorFESpaceFunction;
import distorted.DistortedVectorCellIntegral;
import distorted.DistortedVectorRightHandSideIntegral;
import distorted.DistortedVectorShapeFunction;
import distorted.DistortedVectorSpace;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import tensorproduct.QuadratureRule1D;

import java.util.ArrayList;
import java.util.List;

public class DiskVectorLaplace
{
	public static void main(final String[] args)
	{
		final int polynomialDegree = 2;
		final DistortedVectorCellIntegral gradGrad = new DistortedVectorCellIntegral(1,
		                                                                             DistortedVectorCellIntegral.GRAD_GRAD,
		                                                                             QuadratureRule1D
			                                                                             .fromPolynomialDegree(
				                                                                             polynomialDegree));
		final DistortedVectorRightHandSideIntegral source = new DistortedVectorRightHandSideIntegral(
			LaplaceReferenceSolution.vectorRightHandSide(),
			DistortedVectorRightHandSideIntegral.VALUE);
		//final PlotWindow p = new PlotWindow();
		
		for (int i = 0; i < 4; i++)
		{
			final DistortedVectorSpace circle = new DistortedVectorSpace(new CoordinateVector(2), 1, i);
			
			System.out.println("Cells done");
			circle.assembleCells();
			circle.assembleFunctions(polynomialDegree);
			circle.initializeSystemMatrix();
			circle.initializeRhs();
			System.out.println("System Initialized");
			circle.evaluateCellIntegrals(List.of(gradGrad), List.of(source));
			circle.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
			circle.setBoundaryValues(LaplaceReferenceSolution.vectorReferenceSolution());
			System.out.println(circle.getShapeFunctions().size());
			
			System.out.println("System Filled");
//		final Vector solution = circle.getSystemMatrix().solve(circle.getRhs());
			final IterativeSolver iterativeSolver = new IterativeSolver();
			iterativeSolver.showProgress = false;
			final Vector solution = iterativeSolver.solveGMRES(circle.getSystemMatrix(), circle.getRhs(),
			                                                   1e-9);
			
			final VectorFESpaceFunction<DistortedVectorShapeFunction> solutionFunction =
				new VectorFESpaceFunction<>(circle.getShapeFunctions(), solution);
//			p.addPlot(new ScalarPlot2D(solutionFunction,
//			                           circle.generatePlotPoints(6 * circle.getCells().size()),
//			                           30));
//			p.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.scalarReferenceSolution(),
//			                           circle.generatePlotPoints(6 * circle.getCells().size()), 30));
			System.out.println(ConvergenceOrderEstimator.normL2VecDifference(solutionFunction,
			                                                                 LaplaceReferenceSolution.vectorReferenceSolution(),
			                                                                 circle.generateIsotropicPoints(
				                                                                 50)) / 50 / 50);
		}
	}
}
