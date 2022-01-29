package linalg;

import basic.DoubleCompare;
import org.junit.Test;

import static org.junit.Assert.*;

public class DenseMatrixTest1
{
	
	@Test
	public void testGetters()
	{
		final DenseVector vector;
		final DenseMatrix singularMatrix;
		final DenseMatrix invertibleMatrix;
		vector = DenseVector.vectorFromValues(0, 1, 2, 3);
		singularMatrix = DenseMatrix.squareMatrixFromValues(0, 1, 0, 5);
		invertibleMatrix = DenseMatrix.squareMatrixFromValues(0, 1, 1, 1);
		assertEquals(1, vector.getOrder());
		assertEquals(2, singularMatrix.getOrder());
		assertEquals(2, invertibleMatrix.getOrder());
		assertEquals(4, vector.getLength());
		assertEquals(2, singularMatrix.getCols());
		assertEquals(2, invertibleMatrix.getCols());
		assertEquals(2, singularMatrix.getRows());
		assertEquals(2, invertibleMatrix.getRows());
		assertFalse(vector.isSparse());
		assertFalse(singularMatrix.isSparse());
		assertFalse(invertibleMatrix.isSparse());
		assertEquals(4, vector.size());
		assertEquals(4, singularMatrix.size());
		assertEquals(4, invertibleMatrix.size());
		assertTrue(DoubleCompare.almostEqual(vector.euclidianNorm(), Math.sqrt(1 + 4 + 9)));
	}
	
	@Test
	public void testLarge()
	{
		final DenseMatrix largeDense = new DenseMatrix(50, 50);
		final DenseMatrix largeDense2 = new DenseMatrix(50, 50);
		for (int i = 0; i < 100; i++)
		{
			final int x = (int) (Math.random() * 50);
			final int y = (int) (Math.random() * 50);
			final double v = (Math.random() * 4000);
			largeDense.set(v, y, x);
			largeDense2.set(v, y, x);
		}
		for (int i = 0; i < 50; i++)
		{
			for (int j = 0; j < 50; j++)
				assertTrue(DoubleCompare.almostEqual(largeDense.at(i, j), largeDense2.at(i, j)));
		}
		assertTrue(largeDense.almostEqual(largeDense2));
		for (int i = 0; i < 1000; i++)
		{
			final int x = (int) (Math.random() * 50);
			final int y = (int) (Math.random() * 50);
			final double v = (int) (Math.random() * 4000);
			largeDense.add(v, y, x);
			largeDense2.add(v, y, x);
		}
		assertEquals(largeDense, largeDense2);
		final DenseMatrix largeDense3 = largeDense.add(largeDense);
		assertEquals(largeDense, largeDense2);
		final DenseMatrix largeDense4 = largeDense.mul(0.1);
		final DenseMatrix largeDenseRec =
			largeDense
				.getStrictlyLowerTriangleMatrix()
				.add(largeDense.getDiagonalMatrix())
				.add(largeDense.getStrictlyUpperTriangleMatrix());
		assertFalse(largeDense4.isSparse());
		assertEquals(largeDense, largeDense2);
		for (int i = 0; i < 50; i++)
		{
			for (int j = 0; j < 50; j++)
			{
				assertTrue(DoubleCompare.almostEqual(largeDense3.at(i, j), 2 * largeDense.at(i, j)));
				assertTrue(DoubleCompare.almostEqual(largeDense4.at(i, j), 0.1 * largeDense.at(i, j)));
				assertTrue(DoubleCompare.almostEqual(largeDenseRec.at(i, j), largeDense.at(i, j)));
			}
		}
		final DenseVector largeVector = new DenseVector(50);
		final DenseVector sparseVector = new DenseVector(50);
		for (int i = 0; i < 50; i++)
		{
			final double x = Math.random();
			largeVector.set(x, i);
			sparseVector.set(x, i);
		}
		final DenseVector mul1 = largeDense.mvMul(largeVector);
		final DenseVector mul2 = largeDense3.mvMul(largeVector);
		assertEquals(largeDense, largeDense3.mul(0.5));
		for (int i = 0; i < 50; i++)
		{
			assertTrue(DoubleCompare.almostEqual(mul1.mul(2)
			                                         .at(i), mul2.at(i)));
			assertTrue(DoubleCompare.almostEqual(mul1.at(i) * 2, mul2.at(i)));
		}
		assertTrue(mul1.almostEqual(mul2.mul(0.5)));
		for (int i = 0; i < 50; i++)
		{
			largeVector.add(Math.random(), i);
		}
		assertEquals(largeDense.transpose()
		                       .transpose(), largeDense);
		assertTrue(largeDense.tvMul(largeVector)
		                     .almostEqual(largeDense2.transpose()
		                                             .mvMul(largeVector)
		                                 ));
	}
	
	@Test
	public void testSolve()
	{
		final DenseMatrix largeDense = new DenseMatrix(50, 50);
		final DenseVector largeVector = new DenseVector(50);
		for (int i = 0; i < 100; i++)
		{
			final int x = (int) (Math.random() * 50);
			final int y = (int) (Math.random() * 50);
			final double v = (Math.random() * 4);
			final double w = (Math.random() * 4);
			largeVector.set(w, x);
			largeDense.set(v, y, x);
			largeDense.add(5.2, i % 50, i % 50);
		}
		final DenseMatrix denseInverse = largeDense.inverse();
		assertTrue(denseInverse.mvMul(largeVector)
		                       .almostEqual(largeDense.solve(largeVector)));
		assertTrue(denseInverse.mmMul(largeDense)
		                       .almostEqual(largeDense.mmMul(denseInverse)));
		assertTrue(denseInverse.mmMul(largeDense)
		                       .almostEqual(DenseMatrix.identity(50)));
		assertTrue(denseInverse.inverse()
		                       .almostEqual(largeDense));
		assertTrue(new IterativeSolver()
			           .solveBiCGStab(largeDense, largeVector, 1e-10)
			           .almostEqual(denseInverse.mvMul(largeVector), 1e-8
			                       ));
		assertTrue(largeDense.mvMul(new IterativeSolver().solveBiCGStab(largeDense, largeVector,
		                                                                1e-12))
		                     .almostEqual(largeVector));
		assertTrue(new IterativeSolver()
			           .solveGMRES(largeDense, largeVector, 1e-10)
			           .almostEqual(denseInverse.mvMul(largeVector), 1e-10
			                       ));
		assertTrue(new IterativeSolver()
			           .solveBiCGStab(largeDense, largeVector, 1e-10)
			           .almostEqual(denseInverse.mvMul(largeVector), 1e-10
			                       ));
		assertTrue(new IterativeSolver()
			           .solveGMRES(largeDense, largeVector, 1e-10)
			           .almostEqual(denseInverse.mvMul(largeVector), 1e-10
			                       ));
		final DenseMatrix symmDense = denseInverse.add(denseInverse.transpose());
		assertTrue(symmDense.solveSymm(largeVector)
		                    .almostEqual(new IterativeSolver().solveCG(symmDense,
		                                                               largeVector,
		                                                               1e-10), 1e-10));
		assertTrue(symmDense.solve(largeVector)
		                    .almostEqual(new IterativeSolver().solveCG(symmDense,
		                                                               largeVector, 1e-10),
		                                 1e-10));
		assertTrue(symmDense.solve(largeVector)
		                    .almostEqual(symmDense.solveSymm(largeVector)));
	}
	
	@Test
	public void testSpvecSolve()
	{
		final DenseVector small = new DenseVector(40);
		small.add(4.3, 4);
		small.add(7.8, 5);
		assertTrue(DoubleCompare.almostEqual(small.euclidianNorm(), Math.sqrt(4.3 * 4.3 + 7.8 * 7.8)));
		final DenseMatrix largeDense = new DenseMatrix(50, 50);
		final DenseVector largeVector = new DenseVector(50);
		final DenseVector largeVector2 = new DenseVector(50);
		for (int i = 0; i < 100; i++)
		{
			final int x = (int) (Math.random() * 50);
			final int y = (int) (Math.random() * 50);
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
			largeDense.add(5.2, i % 50, i % 50);
		}
		DenseMatrix denseInverse = largeDense.inverse();
		for (int j = 0; j < 50; j++)
		{
			assertTrue(DoubleCompare.almostEqual(largeVector.at(j), largeVector2.at(j)));
		}
		for (int j = 0; j < 50; j++)
		{
			assertTrue(DoubleCompare.almostEqual(largeVector.mul(0.1)
			                                                .at(j), largeVector2.at(j) * 0.1));
		}
		assertTrue(largeVector.almostEqual(largeVector2));
		assertTrue(largeVector.add(largeVector)
		                      .almostEqual(largeVector2.mul(2)));
		assertTrue(largeVector.sub(largeVector)
		                      .euclidianNorm() < 1e-14);
		assertTrue(largeVector2.sub(largeVector)
		                       .euclidianNorm() < 1e-14);
		assertTrue(largeVector2.sub(largeVector2)
		                       .euclidianNorm() < 1e-14);
		assertTrue(largeVector.sub(largeVector2)
		                      .euclidianNorm() < 1e-14);
		assertTrue(largeDense.mvMul(largeVector)
		                     .almostEqual(largeDense.mvMul(largeVector)));
		assertTrue(largeDense.mvMul(largeVector)
		                     .almostEqual(largeDense.mvMul(largeVector)));
		assertTrue(denseInverse.mvMul(largeVector)
		                       .almostEqual(largeDense.solve(largeVector)));
		assertTrue(denseInverse.mmMul(largeDense)
		                       .almostEqual(largeDense.mmMul(denseInverse)));
		assertTrue(denseInverse.mmMul(largeDense)
		                       .almostEqual(DenseMatrix.identity(50)));
		final SparseMatrix id = SparseMatrix.identity(largeDense.getCols());
		final SparseMatrix sparseInv = new SparseMatrix(denseInverse);
		assertEquals(new SparseMatrix(largeDense), largeDense);
		assertEquals(new DenseMatrix(new SparseMatrix(largeDense)), largeDense);
		assertEquals(new SparseMatrix(denseInverse), denseInverse);
		assertEquals(new DenseMatrix(new SparseMatrix(denseInverse)), denseInverse);
		assertEquals(largeDense.mmMul(id), largeDense);
		assertEquals(sparseInv.mmMul(largeDense), DenseMatrix.identity(50));
		assertEquals(largeDense.mmMul(sparseInv), DenseMatrix.identity(50));
		assertTrue(denseInverse.inverse()
		                       .almostEqual(largeDense));
		denseInverse = largeDense.inverse();
		final DenseVector sparsecopy = new DenseVector(largeVector);
		assertTrue(largeDense.mvMul(new IterativeSolver().solveBiCGStab(largeDense, sparsecopy,
		                                                                1e-12))
		                     .almostEqual(sparsecopy));
		assertTrue(largeDense.mvMul(new IterativeSolver().solveBiCGStab(largeDense, largeVector,
		                                                                1e-12))
		                     .almostEqual(sparsecopy));
		assertTrue(new IterativeSolver()
			           .solveBiCGStab(largeDense, largeVector, 1e-12)
			           .almostEqual(denseInverse.mvMul(largeVector), 1e-10
			                       ));
		assertTrue(new IterativeSolver()
			           .solveGMRES(largeDense, largeVector, 1e-10)
			           .almostEqual(denseInverse.mvMul(largeVector), 1e-10
			                       ));
		assertTrue(new IterativeSolver()
			           .solveBiCGStab(largeDense, largeVector, 1e-10)
			           .almostEqual(denseInverse.mvMul(largeVector), 1e-10
			                       ));
		assertTrue(new IterativeSolver()
			           .solveGMRES(largeDense, largeVector, 1e-10)
			           .almostEqual(denseInverse.mvMul(largeVector), 1e-10
			                       ));
		final DenseMatrix symmDense = denseInverse.add(denseInverse.transpose());
		assertTrue(symmDense.solveSymm(largeVector)
		                    .almostEqual(new IterativeSolver().solveCG(symmDense,
		                                                               largeVector,
		                                                               1e-10), 1e-10));
		assertTrue(symmDense.solve(largeVector)
		                    .almostEqual(new IterativeSolver().solveCG(symmDense,
		                                                               largeVector, 1e-10),
		                                 1e-10));
		assertTrue(symmDense.solve(largeVector)
		                    .almostEqual(symmDense.solveSymm(largeVector)));
	}
}
