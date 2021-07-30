package linalg;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class GMResTest
{
	@Test
	public void testIdentity() {
		SparseMatrix id = SparseMatrix.identity(15);
		DenseVector b = DenseVector.vectorFromValues(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
		IterativeSolver i = new IterativeSolver();
		assertEquals(b,i.solveGMRES(id,b,1e-12));
	}
	@Test
	public void testSymm() {
		DenseMatrix symm = DenseMatrix.squareMatrixFromValues(0,1,2,1,0,3,2,3,0);
		DenseVector b = DenseVector.vectorFromValues(3,4,5);
		IterativeSolver i = new IterativeSolver();
		Vector sol = i.solveGMRES(symm, b, 1e-12);
		assertEquals(b, symm.mvMul(sol));
		assertEquals(symm.solve(b), sol);
	}
	@Test
	public void testNonSymm() {
		DenseMatrix nonsymm = DenseMatrix.squareMatrixFromValues(3,1,7,2,0,3,2,3,0);
		DenseVector b = DenseVector.vectorFromValues(3,4,5);
		IterativeSolver i = new IterativeSolver();
		Vector sol = i.solveGMRES(nonsymm, b, 1e-12);
		assertEquals(b, nonsymm.mvMul(sol));
		assertEquals(nonsymm.solve(b), sol);
	}
	@Test
	public void testLarge() {
		int n = 10000;
		SparseMatrix large = new SparseMatrix(n,n);
		DenseVector b = new DenseVector(n);
		for(int i = 0; i < n*10; i++)
		{
			b.add(1,i%n);
			large.add(0.2,i%n,i%n);
			large.add(Math.random(), (int)(Math.random()*n), (int)(Math.random()*n));
		}
		IterativeSolver i = new IterativeSolver();
		Vector sol = i.solveGMRES(large, b, 1e-10);
		assertTrue(b.almostEqual(large.mvMul(sol)));
	}
	@Test
	public void testPIdentity() {
		SparseMatrix id = SparseMatrix.identity(15);
		DenseVector b = DenseVector.vectorFromValues(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
		IterativeSolver i = new IterativeSolver();
		assertEquals(b,i.solvePGMRES(id, id,b,1e-12));
	}
	@Test
	public void testPSymm() {
		DenseMatrix symm = DenseMatrix.squareMatrixFromValues(0,1,2,1,0,3,2,3,0);
		DenseVector b = DenseVector.vectorFromValues(3,4,5);
		IterativeSolver i = new IterativeSolver();
		Vector sol = i.solvePGMRES(symm, SparseMatrix.identity(3), b, 1e-12);
	}
	@Test
	public void testPNonSymm() {
		DenseMatrix nonsymm = DenseMatrix.squareMatrixFromValues(3,1,7,2,0,3,2,3,0);
		DenseVector b = DenseVector.vectorFromValues(3,4,5);
		IterativeSolver i = new IterativeSolver();
		Vector sol = i.solvePGMRES(nonsymm,SparseMatrix.identity(3), b, 1e-12);
		assertEquals(b, nonsymm.mvMul(sol));
		assertEquals(nonsymm.solve(b), sol);
	}
	@Test
	public void testPLarge() {
		int n = 10000;
		SparseMatrix large = new SparseMatrix(n,n);
		DenseVector b = new DenseVector(n);
		for(int i = 0; i < n*10; i++)
		{
			b.add(1,i%n);
			large.add(0.2,i%n,i%n);
			large.add(Math.random(), (int)(Math.random()*n), (int)(Math.random()*n));
		}
		IterativeSolver i = new IterativeSolver();
		Vector sol = i.solvePGMRES(large, SparseMatrix.identity(n), b, 1e-10);
		assertTrue(b.almostEqual(large.mvMul(sol)));
		
	}
	@Test
	public void testGoodPNonSymm() {
		DenseMatrix nonsymm = DenseMatrix.squareMatrixFromValues(3,1,7,2,4,3,2,3,7);
		DenseVector b = DenseVector.vectorFromValues(3,4,5);
		IterativeSolver i = new IterativeSolver();
		Vector sol = i.solvePGMRES(nonsymm,nonsymm.inverse(), b, 1e-12);
		assertEquals(b, nonsymm.mvMul(sol));
		assertEquals(nonsymm.solve(b), sol);
	}
}
