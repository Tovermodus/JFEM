package tensorproduct;

import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.Matrix;
import org.junit.jupiter.api.Test;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class TestReferenceFaceIntegral
{
	@Test
	public void test2DJumpJumpMatrix()
	{
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -10);
		CoordinateVector end = CoordinateVector.fromValues(4, 1);
		int polynomialDegree = 1;
		
		TPFESpace grid = new TPFESpace(start, end,
			Ints.asList(7, 11));
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
			Ints.asList(7, 11));
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
		
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix()));
		
	}
	@Test
	public void test2DvaluenormalaverageMatrix()
	{
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -10);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 3;
		
		TPFESpace grid = new TPFESpace(start, end,
			Ints.asList(3, 5));
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
			Ints.asList(3, 5));
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
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix()));
		
	}
	
	@Test
	public void test3DJumpJumpMatrix()
	{
		Matrix referenceMatrix;
		CoordinateVector start = CoordinateVector.fromValues(-1, -10,6);
		CoordinateVector end = CoordinateVector.fromValues(1, 1,17);
		int polynomialDegree = 1;
		
		TPFESpace grid = new TPFESpace(start, end,
			Ints.asList(3, 5,2));
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
			Ints.asList(3, 5,2));
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
		assertTrue(referenceMatrix.almostEqual(grid.getSystemMatrix()));
		
	}
	
	@Test
	public void test2DJumpJump()
	{
		CartesianGrid g = new CartesianGrid(CoordinateVector.fromValues(3, 130),
			CoordinateVector.fromValues(6,138),new IntCoordinates(1,2));
		TPFace f = g.faces.get(5);
		assertEquals(f.getNormalUpstreamCell().center(), CoordinateVector.fromValues(4.5, 132.0));
		assertEquals(f.getNormalDownstreamCell().center(), CoordinateVector.fromValues(4.5, 136));
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
