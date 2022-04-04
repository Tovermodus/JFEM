package linalg;

import basic.DoubleCompare;
import com.google.common.base.Stopwatch;
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.*;

public class SparseMatrixTest1
{
	@Test
	public void testLarge()
	{
		final int n = 100;
		final SparseMatrix largeSparse = new SparseMatrix(n, n);
		final DenseMatrix largeDense = new DenseMatrix(n, n);
		for (int i = 0; i < n * 4; i++)
		{
			final int x = (int) (Math.random() * n);
			final int y = (int) (Math.random() * n);
			final double v = (Math.random() * 4000);
			largeSparse.set(v, y, x);
			largeDense.set(v, y, x);
		}
		largeSparse.set(123, 0, 0);
		largeDense.set(123, 0, 0);
		largeSparse.set(1224, 0, 0);
		largeDense.set(1224, 0, 0);
		int nonzeros = 0;
		for (final IntCoordinates c : largeDense.getShape()
		                                        .range())
			if (largeDense.at(c) != 0)
			{
				nonzeros++;
			}
		assertEquals(largeSparse.getCoordinateEntryList()
		                        .size(), nonzeros);
		assertTrue(largeSparse.almostEqual(largeDense));
		for (int i = 0; i < 1000; i++)
		{
			final int x = (int) (Math.random() * n);
			final int y = (int) (Math.random() * n);
			final double v = (int) (Math.random() * 4000);
			largeSparse.add(v, y, x);
			largeDense.add(v, y, x);
			largeSparse.add(v, y, x);
			largeDense.add(v, y, x);
		}
		assertTrue(largeSparse.almostEqual(largeDense));
		assertEquals(largeSparse, largeDense);
		final SparseMatrix largeDense3 = largeSparse.add(largeSparse);
		assertEquals(largeSparse, largeDense);
		final SparseMatrix largeDense4 = largeSparse.mul(0.1);
		final SparseMatrix largeDenseRec = largeSparse
			.getStrictlyLowerTriangleMatrix()
			.add(largeSparse.getDiagonalMatrix())
			.add(largeSparse.getStrictlyUpperTriangleMatrix());
		assertTrue(largeDense4.isSparse());
		assertEquals(largeSparse, largeDense);
		
		assertTrue(largeSparse.almostEqual(largeDense));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				assertTrue(Math.abs(largeDense4.at(i, j) - 0.1 * largeSparse.at(i, j)) < 1e-10);
				assertTrue(Math.abs(largeDense3.at(i, j) - 2 * largeDense.at(i, j)) < 1e-10);
				assertTrue(DoubleCompare.almostEqual(largeDenseRec.at(i, j), largeSparse.at(i, j)));
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
		assertEquals(largeSparse, largeDense3.mul(0.5));
		assertEquals(mul1.mul(2), mul2);
		assertTrue(mul1.almostEqual(mul2.mul(0.5)));
		for (int i = 0; i < n; i++)
		{
			largeVector.add(Math.random(), i);
		}
		assertEquals(largeSparse.transpose()
		                        .transpose(), largeSparse);
		assertTrue(largeSparse.tvMul(largeVector)
		                      .almostEqual(largeDense.transpose()
		                                             .mvMul(largeVector)));
		assertEquals(largeSparse.mmMul(largeSparse), largeSparse.mmMul(largeDense));
		assertEquals(largeSparse.mmMul(largeSparse), largeDense.mmMul(largeSparse));
		assertEquals(largeSparse.mmMul(largeSparse), largeDense.mmMul(largeDense));
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
		assertTrue(denseInverse.mvMul(largeVector)
		                       .almostEqual(largeDense.solve(largeVector)));
		assertTrue(largeDense.mmMul(denseInverse)
		                     .almostEqual(SparseMatrix.identity(50)));
		assertTrue(denseInverse.mmMul(largeDense)
		                       .almostEqual(SparseMatrix.identity(50)));
		assertTrue(denseInverse.mmMul(largeDense)
		                       .almostEqual(largeDense.mmMul(denseInverse)));
		assertTrue(denseInverse.inverse()
		                       .almostEqual(largeDense));
		
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
		
		final SparseMatrix symmDense = largeDense.add(largeDense.transpose());
		assertTrue(symmDense
			           .solveSymm(largeVector)
			           .almostEqual(new IterativeSolver().solveCG(symmDense, largeVector, 1e-11)));
		assertTrue(symmDense
			           .solve(largeVector)
			           .almostEqual(new IterativeSolver().solveCG(symmDense, largeVector, 1e-11)));
		assertTrue(symmDense.solve(largeVector)
		                    .almostEqual(symmDense.solveSymm(largeVector)));
	}
	
	@Test
	public void testEquals()
	{
		final SparseMatrix mat1 = new SparseMatrix(10, 10);
		assertEquals(mat1, new SparseMatrix(10, 10));
		mat1.add(1e-15, 0, 0);
		assertEquals(mat1, new SparseMatrix(10, 10));
		mat1.add(1e-12, 0, 0);
		assertEquals(mat1, new SparseMatrix(10, 10));
		assertEquals(new SparseMatrix(10, 10), mat1);
		mat1.add(1e-9, 0, 0);
		assertNotEquals(mat1, new SparseMatrix(10, 10));
		assertNotEquals(new SparseMatrix(10, 10), mat1);
	}
	
	@Test
	public void testSpvecSolve()
	{
		final DenseVector small = new DenseVector(40);
		small.add(4.3, 4);
		small.add(7.8, 5);
		assertTrue(DoubleCompare.almostEqual(small.euclidianNorm(), Math.sqrt(4.3 * 4.3 + 7.8 * 7.8)));
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
		
		final DenseVector sparsecopy = new DenseVector(largeVector);
		
		assertTrue(largeDense
			           .mvMul(new IterativeSolver().solveBiCGStab(largeDense, sparsecopy, 1e-12))
			           .almostEqual(sparsecopy));
		
		assertTrue(new IterativeSolver()
			           .solveBiCGStab(largeDense, largeVector, 1e-12)
			           .almostEqual(new IterativeSolver().solveGMRES(largeDense, largeVector, 1e-10)));
		final SparseMatrix symmDense = largeDense.add(largeDense.transpose());
		
		assertTrue(symmDense
			           .solveSymm(largeVector)
			           .almostEqual(new IterativeSolver().solveCG(symmDense, largeVector, 1e-10)));
	}
	
	@Test
	public void testmmMul()
	{
		final int n = 2002;
		final DenseMatrix dense = new DenseMatrix(n, n);
		final SparseMatrix sparse = new SparseMatrix(n, n);
		final Random generator = new Random(3145);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 50; j++)
				sparse.add(generator.nextDouble(), i, generator.nextInt(n));
			for (int j = 0; j < n; j++)
				dense.add(generator.nextDouble(), i, j);
		}
		final DenseMatrix sparseToDense = new DenseMatrix(sparse);
		final Stopwatch s = Stopwatch.createStarted();
		DenseMatrix result = null;
		for (int i = 0; i < 10; i++)
		{
			result = dense.mmMul(sparse);
		}
		System.out.println(s.elapsed());
		assertEquals(result, dense.mmMul(sparseToDense));
		System.out.println(s.elapsed());
		result = sparseToDense.mmMul(sparse);
		System.out.println(s.elapsed());
		final SparseMatrix spspresult = sparse.mmMul(sparse);
		System.out.println(s.elapsed());
		assertEquals(spspresult, result);
		System.out.println(s.elapsed());
		assertEquals(spspresult, sparse.mmMul(sparseToDense));
		System.out.println(s.elapsed());
	}
	
	public static SparseMatrix createLargeSparseMatrix(final int n)
	{
		
		final Random generator = new Random(1627896);
		final SparseMatrix ret = new SparseMatrix(n, n);
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			if (c.get(0) == c.get(1))
				ret.add(100, c);
			if (generator.nextDouble() < 0.05)
				ret.add(generator.nextDouble() - 0.4, c);
		}
		return ret;
	}
	
	public static DenseVector createRhs(final Matrix m)
	{
		
		final Random generator = new Random(323781645);
		final DenseVector ret = new DenseVector(m.getVectorSize());
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			ret.add(generator.nextDouble() - 0.4, c);
		}
		return ret;
	}
	
	@Test
	public void testLLarge()
	{
		SparseMatrix s = createLargeSparseMatrix(2000);
		final DenseVector b = createRhs(s);
		System.out.println("created");
		Stopwatch st = Stopwatch.createStarted();
		final Vector x1 = s.solve(b);
		System.out.println(st.elapsed());
		st = Stopwatch.createStarted();
		final Vector x2 = new DenseVector(s
			                                  .getSparseSolver()
			                                  .solve(b.entries), true);
		System.out.println(st.elapsed());
		assertEquals(x1, x2);
		s = createLargeSparseMatrix(2000);
		st = Stopwatch.createStarted();
		final DenseMatrix d1 = s.inverse();
		System.out.println(st.elapsed());
		st = Stopwatch.createStarted();
		final DenseMatrix d2 = new DenseMatrix(s.getSparseSolver()
		                                        .inverse());
		System.out.println(st.elapsed());
		assertEquals(d1, d2);
	}
}
