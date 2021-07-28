package tensorproduct;

import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.Matrix;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class TestReferenceFaceIntegral
{
	@Test
	public void test2DJumpJumpMatrix()
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		builder.build();
		long startTime = System.nanoTime();
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -10);
		CoordinateVector end = CoordinateVector.fromValues(4, 1);
		int polynomialDegree = 1;
		
		TPFESpace grid = new TPFESpace(start, end,
			Ints.asList(7, 11), polynomialDegree);
		TPFaceIntegral<TPShapeFunction> jj = new TPFaceIntegral<>(ScalarFunction.constantFunction(1),
			TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
		ArrayList<FaceIntegral<TPFace, TPShapeFunction>> faceIntegrals =
			new ArrayList<>();
		faceIntegrals.add(jj);
		ArrayList<BoundaryRightHandSideIntegral<TPFace, TPShapeFunction>> rightHandSideIntegrals =
			new ArrayList<>();
		grid.assembleCells();
		System.out.println(grid.getCells());
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateFaceIntegrals(faceIntegrals, rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new TPFESpace(start, end,
			Ints.asList(7, 11), polynomialDegree);
		jj = new TPFaceIntegralViaReferenceFace<>(1,
			TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
		faceIntegrals =
			new ArrayList<>();
		faceIntegrals.add(jj);
		rightHandSideIntegrals = new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateFaceIntegrals(faceIntegrals, rightHandSideIntegrals);
		
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix(), 1e-12));
		
	}
	@Test
	public void test2DvaluenormalaverageMatrix()
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		long startTime = System.nanoTime();
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -10);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 3;
		
		TPFESpace grid = new TPFESpace(start, end,
			Ints.asList(3, 5), polynomialDegree);
		TPFaceIntegral<TPShapeFunction> jj = new TPFaceIntegral<>(ScalarFunction.constantFunction(1),
			TPFaceIntegral.VALUE_JUMP_GRAD_NORMALAVERAGE);
		ArrayList<FaceIntegral<TPFace, TPShapeFunction>> faceIntegrals =
			new ArrayList<>();
		faceIntegrals.add(jj);
		ArrayList<BoundaryRightHandSideIntegral<TPFace, TPShapeFunction>> rightHandSideIntegrals =
			new ArrayList<>();
		grid.assembleCells();
		System.out.println(grid.getCells());
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateFaceIntegrals(faceIntegrals, rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new TPFESpace(start, end,
			Ints.asList(3, 5), polynomialDegree);
		jj = new TPFaceIntegralViaReferenceFace<>(1,
			TPFaceIntegral.VALUE_JUMP_GRAD_NORMALAVERAGE);
		faceIntegrals =
			new ArrayList<>();
		faceIntegrals.add(jj);
		rightHandSideIntegrals = new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateFaceIntegrals(faceIntegrals, rightHandSideIntegrals);
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix(), 1e-12));
		
	}
	
	@Test
	public void test3DJumpJumpMatrix()
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		long startTime = System.nanoTime();
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -10,6);
		CoordinateVector end = CoordinateVector.fromValues(1, 1,17);
		int polynomialDegree = 1;
		
		TPFESpace grid = new TPFESpace(start, end,
			Ints.asList(3, 5,2), polynomialDegree);
		TPFaceIntegral<TPShapeFunction> jj = new TPFaceIntegral<>(ScalarFunction.constantFunction(1),
			TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
		ArrayList<FaceIntegral<TPFace, TPShapeFunction>> faceIntegrals =
			new ArrayList<>();
		faceIntegrals.add(jj);
		ArrayList<BoundaryRightHandSideIntegral<TPFace, TPShapeFunction>> rightHandSideIntegrals =
			new ArrayList<>();
		grid.assembleCells();
		System.out.println(grid.getCells());
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateFaceIntegrals(faceIntegrals, rightHandSideIntegrals);
		referenceMatrix = grid.getSystemMatrix();
		
		grid = new TPFESpace(start, end,
			Ints.asList(3, 5,2), polynomialDegree);
		jj = new TPFaceIntegralViaReferenceFace<>(1,
			TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
		faceIntegrals =
			new ArrayList<>();
		faceIntegrals.add(jj);
		rightHandSideIntegrals = new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		grid.evaluateFaceIntegrals(faceIntegrals, rightHandSideIntegrals);
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix(), 1e-12));
		
	}
	
	@Test
	public void test2DJumpJump()
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.parallelizeThreads = false;
		ArrayList<Cell1D> cell1DList = new ArrayList<>();
		cell1DList.add(new Cell1D(3, 6));
		TPFace f = new TPFace(cell1DList,1,134, false);
		cell1DList.add(new Cell1D(130,134));
		f.addCell(new TPCell(cell1DList));
		cell1DList.set(1,new Cell1D(134,137));
		f.addCell(new TPCell(cell1DList));
		assertEquals(f.getNormalUpstreamCell().center(), CoordinateVector.fromValues(4.5, 132.0));
		assertEquals(f.getNormalDownstreamCell().center(), CoordinateVector.fromValues(4.5, 135.5));
		for(int k = 0; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					TPShapeFunction f1 = new TPShapeFunction(f.getNormalDownstreamCell(), k, i);
					TPShapeFunction f2 = new TPShapeFunction(f.getNormalUpstreamCell(), k, j);
					TPFaceIntegral<TPShapeFunction> gg = new TPFaceIntegral<>(ScalarFunction.constantFunction(1),
						TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
					TPFaceIntegral<TPShapeFunction> gg2 = new TPFaceIntegralViaReferenceFace<>(1,
						TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
					assertTrue(Math.abs(gg.evaluateFaceIntegral(f, f1, f2) - gg2.evaluateFaceIntegral(f,
						f1, f2)) < 1e-12, "Difference of "+Math.abs(gg.evaluateFaceIntegral(f
						, f1,
						f2) - gg2.evaluateFaceIntegral(f,
						f1, f2))+ " for polynomial degree " + k + " and functions " + i+" and " + j);
				}
		for(int k = 0; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					TPShapeFunction f1 = new TPShapeFunction(f.getNormalUpstreamCell(), k, i);
					TPShapeFunction f2 = new TPShapeFunction(f.getNormalUpstreamCell(), k, j);
					TPFaceIntegral<TPShapeFunction> gg = new TPFaceIntegral<>(ScalarFunction.constantFunction(1),
						TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
					TPFaceIntegral<TPShapeFunction> gg2 = new TPFaceIntegralViaReferenceFace<>(1,
						TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
					assertTrue(Math.abs(gg.evaluateFaceIntegral(f, f1, f2) - gg2.evaluateFaceIntegral(f,
						f1, f2)) < 1e-12, "Difference of "+Math.abs(gg.evaluateFaceIntegral(f
						, f1,
						f2) - gg2.evaluateFaceIntegral(f,
						f1, f2))+ " for polynomial degree " + k + " and functions " + i+" and " + j);
				}
		for(int k = 0; k < 5; k++)
			for (int i = 0; i < (k + 1) * (k + 1); i++)
				for (int j = 0; j < (k + 1) * (k + 1); j++)
				{
					TPShapeFunction f1 = new TPShapeFunction(f.getNormalDownstreamCell(), k, i);
					TPShapeFunction f2 = new TPShapeFunction(f.getNormalUpstreamCell(), k, j);
					TPFaceIntegral<TPShapeFunction> gg = new TPFaceIntegral<>(ScalarFunction.constantFunction(1),
						TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
					TPFaceIntegral<TPShapeFunction> gg2 = new TPFaceIntegralViaReferenceFace<>(1,
						TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
					assertTrue(Math.abs(gg.evaluateFaceIntegral(f, f1, f2) - gg2.evaluateFaceIntegral(f,
						f1, f2)) < 1e-12, "Difference of "+Math.abs(gg.evaluateFaceIntegral(f
						, f1,
						f2) - gg2.evaluateFaceIntegral(f,
						f1, f2))+ " for polynomial degree " + k + " and functions " + i+" and " + j);
				}
	}
}
