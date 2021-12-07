package examples;

import basic.VectorFESpaceFunction;
import distorted.*;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import org.junit.Test;
import tensorproduct.QuadratureRule1D;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class VectorDiskLaplace
{
	@Test(timeout = 30000)
	public void testConvergence()
	{
		final int polynomialDegree = 3;
		final DistortedVectorCellIntegral gradGrad = new DistortedVectorCellIntegral(1,
		                                                                             DistortedVectorCellIntegral.GRAD_GRAD,
		                                                                             QuadratureRule1D
			                                                                             .fromPolynomialDegree(
				                                                                             polynomialDegree));
		final DistortedVectorRightHandSideIntegral source = new DistortedVectorRightHandSideIntegral(
			LaplaceReferenceSolution.vectorRightHandSide(),
			DistortedVectorRightHandSideIntegral.VALUE);
		final CircleVectorSpace circle = new CircleVectorSpace(new CoordinateVector(2), 1, 3);
		
		System.out.println("Cells done");
		circle.assembleCells();
		circle.assembleFunctions(polynomialDegree);
		circle.initializeSystemMatrix();
		circle.initializeRhs();
		System.out.println("System Initialized");
		circle.evaluateCellIntegrals(List.of(gradGrad), List.of(source));
		circle.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		circle.setBoundaryValues(LaplaceReferenceSolution.vectorReferenceSolution());
		System.out.println(circle.getShapeFunctions()
		                         .size());
		
		System.out.println("System Filled");
//		final Vector solution = circle.getSystemMatrix().solve(circle.getRhs());
		final IterativeSolver iterativeSolver = new IterativeSolver();
		iterativeSolver.showProgress = true;
		final Vector solution = iterativeSolver.solveCG(circle.getSystemMatrix(), circle.getRhs(),
		                                                1e-10);
		
		final VectorFESpaceFunction<DistortedVectorShapeFunction> solutionFunction =
			new VectorFESpaceFunction<>(circle.getShapeFunctions(), solution);
		final double norm = ConvergenceOrderEstimator
			.normL2VecDifference(solutionFunction,
			                     LaplaceReferenceSolution.vectorReferenceSolution(),
			                     circle.generatePlotPoints(6 *
				                                               circle
					                                               .getCells()
					                                               .size()));
		System.out.println(norm);
		assertTrue("" + norm, norm < 0.1);
	}
	
	@Test(timeout = 30000)
	public void testRingConvergence()
	{
		final int polynomialDegree = 1;
		final DistortedVectorCellIntegral gradGrad = new DistortedVectorCellIntegral(1,
		                                                                             DistortedVectorCellIntegral.GRAD_GRAD,
		                                                                             QuadratureRule1D
			                                                                             .fromPolynomialDegree(
				                                                                             polynomialDegree));
		final DistortedVectorRightHandSideIntegral source = new DistortedVectorRightHandSideIntegral(
			LaplaceReferenceSolution.vectorRightHandSide(),
			DistortedVectorRightHandSideIntegral.VALUE);
		final RingVectorSpace circle = new RingVectorSpace(new CoordinateVector(2), 0.6, 0.9, 3);
		
		System.out.println("Cells done");
		circle.assembleCells();
		circle.assembleFunctions(polynomialDegree);
		circle.initializeSystemMatrix();
		circle.initializeRhs();
		System.out.println("System Initialized");
		circle.evaluateCellIntegrals(List.of(gradGrad), List.of(source));
		circle.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		circle.setBoundaryValues(LaplaceReferenceSolution.vectorReferenceSolution());
		System.out.println(circle.getShapeFunctions()
		                         .size());
		
		System.out.println("System Filled");
//		final Vector solution = circle.getSystemMatrix().solve(circle.getRhs());
		final IterativeSolver iterativeSolver = new IterativeSolver();
		iterativeSolver.showProgress = true;
		final Vector solution = iterativeSolver.solveCG(circle.getSystemMatrix(), circle.getRhs(),
		                                                1e-10);
		
		final VectorFESpaceFunction<DistortedVectorShapeFunction> solutionFunction =
			new VectorFESpaceFunction<>(circle.getShapeFunctions(), solution);
		final double norm = ConvergenceOrderEstimator
			.normL2VecDifference(solutionFunction,
			                     LaplaceReferenceSolution.vectorReferenceSolution(),
			                     circle.generateIsotropicPoints(100));
		System.out.println(norm);
		assertTrue("" + norm, norm < 0.05);
	}
}
