package linalg;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class BiCGStabTest
{
	@Test
	public void testIdentity()
	{
		final SparseMatrix id = SparseMatrix.identity(5);
		final DenseVector b = DenseVector.vectorFromValues(1, 2, 3, 4, 5);
		final IterativeSolver i = new IterativeSolver();
		assertEquals(b, i.solveBiCGStab(id, b, 1e-12));
	}
	
	@Test
	public void testSymm()
	{
		final DenseMatrix symm = DenseMatrix.squareMatrixFromValues(0, 1, 2, 1, 0, 3, 2, 3, 0);
		final DenseVector b = DenseVector.vectorFromValues(3, 4, 5);
		final IterativeSolver i = new IterativeSolver();
		final Vector sol = i.solveBiCGStab(symm, b, 1e-12);
		assertEquals(b, symm.mvMul(sol));
		assertEquals(symm.solve(b), sol);
	}
	
	@Test
	public void testNonSymm()
	{
		final DenseMatrix nonsymm = DenseMatrix.squareMatrixFromValues(3, 1, 7, 2, 0, 3, 2, 3, 0);
		final DenseVector b = DenseVector.vectorFromValues(3, 4, 5);
		final IterativeSolver i = new IterativeSolver();
		final Vector sol = i.solveBiCGStab(nonsymm, b, 1e-12);
		assertEquals(b, nonsymm.mvMul(sol));
		assertEquals(nonsymm.solve(b), sol);
	}
	
	@Test
	public void testLarge()
	{
		final int n = 10000;
		final SparseMatrix large = new SparseMatrix(n, n);
		final DenseVector b = new DenseVector(n);
		for (int i = 0; i < n * 10; i++)
		{
			b.add(1, i % n);
			large.add(0.2, i % n, i % n);
			large.add(Math.random(), (int) (Math.random() * n), (int) (Math.random() * n));
		}
		final IterativeSolver i = new IterativeSolver();
		final Vector sol = i.solveBiCGStab(large, b, 1e-10);
		assertTrue(b.almostEqual(large.mvMul(sol)));
	}
	
	@Test
	public void testPIdentity()
	{
		final SparseMatrix id = SparseMatrix.identity(15);
		final DenseVector b = DenseVector.vectorFromValues(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
		final IterativeSolver i = new IterativeSolver();
		assertEquals(b, i.solvePBiCGStab(id, id, b, 1e-12));
	}
	
	@Test
	public void testPSymm()
	{
		final DenseMatrix symm = DenseMatrix.squareMatrixFromValues(0, 1, 2, 1, 0, 3, 2, 3, 0);
		final DenseVector b = DenseVector.vectorFromValues(3, 4, 5);
		final IterativeSolver i = new IterativeSolver();
		final Vector sol = i.solvePBiCGStab(symm, SparseMatrix.identity(3), b, 1e-12);
	}
	
	@Test
	public void testPNonSymm()
	{
		final DenseMatrix nonsymm = DenseMatrix.squareMatrixFromValues(3, 1, 7, 2, 0, 3, 2, 3, 0);
		final DenseVector b = DenseVector.vectorFromValues(3, 4, 5);
		final IterativeSolver i = new IterativeSolver();
		final Vector sol = i.solvePBiCGStab(nonsymm, SparseMatrix.identity(3), b, 1e-12);
		assertEquals(b, nonsymm.mvMul(sol));
		assertEquals(nonsymm.solve(b), sol);
	}
	
	@Test
	public void testPPNonSymm()
	{
		final DenseMatrix nonsymm = DenseMatrix.squareMatrixFromValues(3, 1, 7, 2, 0, 3, 2, 3, 0);
		final DenseVector b = DenseVector.vectorFromValues(3, 4, 5);
		final IterativeSolver i = new IterativeSolver();
		final DenseMatrix precond = nonsymm.inverse();
		precond.add(1, 0, 0);
		final Vector sol = i.solvePBiCGStab(nonsymm, precond, b, 1e-12);
		assertEquals(b, nonsymm.mvMul(sol));
		assertEquals(nonsymm.solve(b), sol);
	}
	
	@Test
	public void testPLarge()
	{
		final int n = 10000;
		final SparseMatrix large = new SparseMatrix(n, n);
		final DenseVector b = new DenseVector(n);
		for (int i = 0; i < n * 10; i++)
		{
			b.add(1, i % n);
			large.add(0.2, i % n, i % n);
			large.add(Math.random(), (int) (Math.random() * n), (int) (Math.random() * n));
		}
		final IterativeSolver i = new IterativeSolver();
		final Vector sol = i.solvePBiCGStab(large, SparseMatrix.identity(n), b, 1e-10);
		assertTrue(b.almostEqual(large.mvMul(sol)));
	}
	
	@Test
	public void testGoodPNonSymm()
	{
		final DenseMatrix nonsymm = DenseMatrix.squareMatrixFromValues(3, 1, 7, 2, 4, 3, 2, 3, 7);
		final DenseVector b = DenseVector.vectorFromValues(3, 4, 5);
		final IterativeSolver i = new IterativeSolver();
		final Vector sol = i.solvePBiCGStab(nonsymm, nonsymm.inverse(), b, 1e-12);
		assertEquals(b, nonsymm.mvMul(sol));
		assertEquals(nonsymm.solve(b), sol);
	}
}
