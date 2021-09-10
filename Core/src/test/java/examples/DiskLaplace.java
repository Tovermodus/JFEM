package examples;

import basic.ScalarFESpaceFunction;
import distorted.DistortedCellIntegral;
import distorted.DistortedRightHandSideIntegral;
import distorted.DistortedShapeFunction;
import distorted.DistortedSpace;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import org.junit.Test;
import tensorproduct.QuadratureRule1D;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class DiskLaplace
{
	@Test
	public void testConvergence()
	{
		final int polynomialDegree = 3;
		final DistortedCellIntegral gradGrad = new DistortedCellIntegral(1, DistortedCellIntegral.GRAD_GRAD,
		                                                                 QuadratureRule1D
			                                                                 .fromPolynomialDegree(
				                                                                 polynomialDegree));
		final DistortedRightHandSideIntegral source = new DistortedRightHandSideIntegral(
			LaplaceReferenceSolution.scalarRightHandSide(),
			DistortedRightHandSideIntegral.VALUE);
		final DistortedSpace circle = new DistortedSpace(new CoordinateVector(2), 1, 3);
		
		System.out.println("Cells done");
		circle.assembleCells();
		circle.assembleFunctions(polynomialDegree);
		circle.initializeSystemMatrix();
		circle.initializeRhs();
		System.out.println("System Initialized");
		circle.evaluateCellIntegrals(List.of(gradGrad), List.of(source));
		circle.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		circle.setBoundaryValues(LaplaceReferenceSolution.scalarReferenceSolution());
		System.out.println(circle.getShapeFunctions().size());
		
		System.out.println("System Filled");
//		final Vector solution = circle.getSystemMatrix().solve(circle.getRhs());
		final IterativeSolver iterativeSolver = new IterativeSolver();
		iterativeSolver.showProgress = false;
		final Vector solution = iterativeSolver.solveGMRES(circle.getSystemMatrix(), circle.getRhs(),
		                                                   1e-10);
		
		final ScalarFESpaceFunction<DistortedShapeFunction> solutionFunction =
			new ScalarFESpaceFunction<>(circle.getShapeFunctions(), solution);
		assertTrue(ConvergenceOrderEstimator
			           .normL2Difference(solutionFunction,
			                             LaplaceReferenceSolution.scalarReferenceSolution(),
			                             circle.generatePlotPoints(6 *
				                                                       circle
					                                                       .getCells()
					                                                       .size()))
			           < 0.05);
	}
}
