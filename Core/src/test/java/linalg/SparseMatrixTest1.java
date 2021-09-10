package linalg;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class SparseMatrixTest1
{
	@Test
	public void testLarge()
	{
		final int n = 100;
		final SparseMatrix largeSparse = new SparseMatrix(n, n);
		final DenseMatrix largeDense2 = new DenseMatrix(n, n);
		for (int i = 0; i < n * 4; i++)
		{
			final int x = (int) (Math.random() * n);
			final int y = (int) (Math.random() * n);
			final double v = (Math.random() * 4000);
			largeSparse.set(v, y, x);
			largeDense2.set(v, y, x);
		}
		largeSparse.set(123, 0, 0);
		largeDense2.set(123, 0, 0);
		largeSparse.set(1224, 0, 0);
		largeDense2.set(1224, 0, 0);
		int nonzeros = 0;
		for (final IntCoordinates c : largeDense2.getShape().range())
			if (largeDense2.at(c) != 0)
			{
				nonzeros++;
			}
		assertEquals(largeSparse.getCoordinateEntryList().size(), nonzeros);
		assertTrue(largeSparse.almostEqual(largeDense2));
		for (int i = 0; i < 1000; i++)
		{
			final int x = (int) (Math.random() * n);
			final int y = (int) (Math.random() * n);
			final double v = (int) (Math.random() * 4000);
			largeSparse.add(v, y, x);
			largeDense2.add(v, y, x);
			largeSparse.add(v, y, x);
			largeDense2.add(v, y, x);
		}
		assertTrue(largeSparse.almostEqual(largeDense2));
		assertEquals(largeSparse, largeDense2);
		final SparseMatrix largeDense3 = (SparseMatrix) largeSparse.add(largeSparse);
		assertEquals(largeSparse, largeDense2);
		final SparseMatrix largeDense4 = largeSparse.mul(0.1);
		final SparseMatrix largeDenseRec = (SparseMatrix) largeSparse
			.getLowerTriangleMatrix()
			.add(largeSparse.getDiagonalMatrix())
			.add(largeSparse.getUpperTriangleMatrix());
		assertTrue(largeDense4.isSparse());
		assertEquals(largeSparse, largeDense2);
		
		assertTrue(largeSparse.almostEqual(largeDense2));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				assertTrue(Math.abs(largeDense4.at(i, j) - 0.1 * largeSparse.at(i, j)) < 1e-10);
				assertTrue(Math.abs(largeDense3.at(i, j) - 2 * largeDense2.at(i, j)) < 1e-10);
				assertEquals(largeDenseRec.at(i, j), largeSparse.at(i, j));
			}
		}
		final DenseVector largeVector = new DenseVector(n);
		final DenseVector DenseVector = new DenseVector(n);
		for (int i = 0; i < n; i++)
		{
			final double x = Math.random();
			largeVector.set(x, i);
			DenseVector.set(x, i);
		}
		final Vector mul1 = largeSparse.mvMul(largeVector);
		final Vector mul2 = largeDense3.mvMul(largeVector);
		final Vector smul1 = largeSparse.mvMul(DenseVector);
		final Vector smul2 = largeDense3.mvMul(DenseVector);
		assertEquals(largeSparse, largeDense3.mul(0.5));
		assertEquals(mul1.mul(2), mul2);
		assertTrue(mul1.almostEqual(mul2.mul(0.5)));
		for (int i = 0; i < n; i++)
		{
			largeVector.add(Math.random(), i);
		}
		assertEquals(largeSparse.transpose().transpose(), largeSparse);
		assertTrue(largeSparse.tvMul(largeVector).almostEqual(largeDense2.transpose().mvMul(largeVector)));
		final Matrix m1 = largeSparse.mmMul(largeDense2);
		final Matrix m2 = largeDense2.mmMul(largeSparse);
		assertTrue(m1.almostEqual(m2));
	}
	
	@Test
	public void testSolve()
	{
		final SparseMatrix largeDense = new SparseMatrix(50, 50);
		final DenseVector largeVector = new DenseVector(50);
		for (int i = 0; i < 4000; i++)
		{
			final int x = (int) (Math.random() * 50);
			final int y = (int) (Math.random() * 50);
			final double v = (Math.random() * 4);
			final double w = (Math.random() * 4);
			largeVector.set(w, x);
			largeDense.set(v, y, x);
			largeDense.add(5.2, i % 50, i % 50);
		}
		
		final DirectlySolvable denseInverse = largeDense.inverse();
		assertTrue(denseInverse.mvMul(largeVector).almostEqual(largeDense.solve(largeVector)));
		assertTrue(largeDense.mmMul(denseInverse).almostEqual(SparseMatrix.identity(50)));
		assertTrue(denseInverse.mmMul(largeDense).almostEqual(SparseMatrix.identity(50)));
		assertTrue(denseInverse.mmMul(largeDense).almostEqual(largeDense.mmMul(denseInverse)));
		assertTrue(denseInverse.inverse().almostEqual(largeDense));
		
		assertTrue(new IterativeSolver()
			           .solveBiCGStab(largeDense, largeVector, 1e-12)
			           .almostEqual(denseInverse.mvMul(largeVector), 1e-10));
		assertTrue(largeDense
			           .mvMul(new IterativeSolver().solveBiCGStab(largeDense, largeVector, 1e-12))
			           .almostEqual(largeVector));
		assertTrue(new IterativeSolver()
			           .solveGMRES(largeDense, largeVector, 1e-12)
			           .almostEqual(denseInverse.mvMul(largeVector), 1e-12));
		assertTrue(new IterativeSolver()
			           .solveBiCGStab(largeDense, largeVector, 1e-12)
			           .almostEqual(denseInverse.mvMul(largeVector)));
		assertTrue(new IterativeSolver()
			           .solveGMRES(largeDense, largeVector, 1e-12)
			           .almostEqual(denseInverse.mvMul(largeVector)));
		
		final SparseMatrix symmDense = (SparseMatrix) largeDense.add(largeDense.transpose());
		assertTrue(symmDense
			           .solveSymm(largeVector)
			           .almostEqual(new IterativeSolver().solveCG(symmDense, largeVector, 1e-12)));
		assertTrue(symmDense
			           .solve(largeVector)
			           .almostEqual(new IterativeSolver().solveCG(symmDense, largeVector, 1e-12)));
		assertTrue(symmDense.solve(largeVector).almostEqual(symmDense.solveSymm(largeVector)));
	}
	
	@Test
	public void testSpvecSolve()
	{
		final DenseVector small = new DenseVector(40);
		small.add(4.3, 4);
		small.add(7.8, 5);
		assertEquals(small.euclidianNorm(), Math.sqrt(4.3 * 4.3 + 7.8 * 7.8));
		final SparseMatrix largeDense = new SparseMatrix(500, 500);
		final DenseVector largeVector = new DenseVector(500);
		final DenseVector largeVector2 = new DenseVector(500);
		for (int i = 0; i < 4000; i++)
		{
			final int x = (int) (Math.random() * 500);
			final int y = (int) (Math.random() * 500);
			final double v = (Math.random() * 4);
			if (i % 10 == 0)
			{
				final double w = (Math.random() * 4);
				largeVector.set(w, x);
				largeVector2.set(w, x);
			}
			if (i % 15 == 0)
			{
				final double w = (Math.random() * 4);
				largeVector.add(w, x);
				largeVector2.add(w, x);
			}
			largeDense.set(v, y, x);
			largeDense.add(5.2, i % 500, i % 500);
		}
		
		assertTrue(largeVector.almostEqual(largeVector2));
		assertTrue(largeVector.add(largeVector).almostEqual(largeVector2.mul(2)));
		assertTrue(largeVector.sub(largeVector).euclidianNorm() < 1e-14);
		assertTrue(largeVector2.sub(largeVector).euclidianNorm() < 1e-14);
		assertTrue(largeVector2.sub(largeVector2).euclidianNorm() < 1e-14);
		assertTrue(largeVector.sub(largeVector2).euclidianNorm() < 1e-14);
		assertTrue(largeDense.mvMul(largeVector).almostEqual(largeDense.mvMul(largeVector)));
		assertTrue(largeDense.mvMul(largeVector).almostEqual(largeDense.mvMul(largeVector)));
		
		final DenseVector sparsecopy = new DenseVector(largeVector);
		
		assertTrue(largeDense
			           .mvMul(new IterativeSolver().solveBiCGStab(largeDense, sparsecopy, 1e-12))
			           .almostEqual(sparsecopy));
		
		assertTrue(new IterativeSolver()
			           .solveBiCGStab(largeDense, largeVector, 1e-12)
			           .almostEqual(new IterativeSolver().solveGMRES(largeDense, largeVector, 1e-10)));
		final SparseMatrix symmDense = (SparseMatrix) largeDense.add(largeDense.transpose());
		
		assertTrue(symmDense
			           .solveSymm(largeVector)
			           .almostEqual(new IterativeSolver().solveCG(symmDense, largeVector, 1e-10)));
	}
}
