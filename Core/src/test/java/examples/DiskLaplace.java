package examples;

import basic.ScalarFESpaceFunction;
import distorted.*;
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
	@Test(timeout = 30000)
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
		final CircleSpace circle = new CircleSpace(new CoordinateVector(2), 1, 3);
		
		System.out.println("Cells done");
		circle.assembleCells();
		circle.assembleFunctions(polynomialDegree);
		circle.initializeSystemMatrix();
		circle.initializeRhs();
		System.out.println("System Initialized");
		circle.evaluateCellIntegrals(List.of(gradGrad), List.of(source));
		circle.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		circle.setBoundaryValues(LaplaceReferenceSolution.scalarReferenceSolution());
		System.out.println(circle.getShapeFunctionMap()
		                         .size());
		
		System.out.println("System Filled");
//		final Vector solution = circle.getSystemMatrix().solve(circle.getRhs());
		final IterativeSolver iterativeSolver = new IterativeSolver();
		iterativeSolver.showProgress = true;
		final Vector solution = iterativeSolver.solveCG(circle.getSystemMatrix(), circle.getRhs(),
		                                                1e-10);
		
		final ScalarFESpaceFunction<DistortedShapeFunction> solutionFunction =
			new ScalarFESpaceFunction<>(circle.getShapeFunctionMap(), solution);
		final double norm = ConvergenceOrderEstimator
			.normL2Difference(solutionFunction,
			                  LaplaceReferenceSolution.scalarReferenceSolution(),
			                  circle.generatePlotPoints(6 *
				                                            circle
					                                            .getCells()
					                                            .size()));
		System.out.println(norm);
		assertTrue("" + norm, norm < 0.05);
	}
	
	@Test(timeout = 100000)
	public void testRingConvergence()
	{
		final int polynomialDegree = 1;
		final DistortedCellIntegral gradGrad = new DistortedCellIntegral(1, DistortedCellIntegral.GRAD_GRAD,
		                                                                 QuadratureRule1D
			                                                                 .fromPolynomialDegree(
				                                                                 polynomialDegree));
		final DistortedRightHandSideIntegral source = new DistortedRightHandSideIntegral(
			LaplaceReferenceSolution.scalarRightHandSide(),
			DistortedRightHandSideIntegral.VALUE);
		final RingSpace circle = new RingSpace(new CoordinateVector(2), 0.6, 0.9, 3);
		
		System.out.println("Cells done");
		circle.assembleCells();
		circle.assembleFunctions(polynomialDegree);
		circle.initializeSystemMatrix();
		circle.initializeRhs();
		System.out.println("System Initialized");
		circle.evaluateCellIntegrals(List.of(gradGrad), List.of(source));
		circle.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		circle.setBoundaryValues(LaplaceReferenceSolution.scalarReferenceSolution());
		System.out.println(circle.getShapeFunctionMap()
		                         .size());
		
		System.out.println("System Filled");
//		final Vector solution = circle.getSystemMatrix().solve(circle.getRhs());
		final IterativeSolver iterativeSolver = new IterativeSolver();
		iterativeSolver.showProgress = true;
		final Vector solution = iterativeSolver.solveCG(circle.getSystemMatrix(), circle.getRhs(),
		                                                1e-10);
		
		final ScalarFESpaceFunction<DistortedShapeFunction> solutionFunction =
			new ScalarFESpaceFunction<>(circle.getShapeFunctionMap(), solution);
		final double norm = ConvergenceOrderEstimator
			.normL2Difference(solutionFunction,
			                  LaplaceReferenceSolution.scalarReferenceSolution(),
			                  circle.generateIsotropicPoints(100));
		System.out.println(norm);
		assertTrue("" + norm, norm < 0.01);
	}
}
