package linalg;

import com.google.common.primitives.Doubles;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;
public class SparseMatrixTest1
{
	@Test
	public void testLarge()
	{
		SparseMatrix largeDense = new SparseMatrix(100,100);
		DenseMatrix largeDense2 = new DenseMatrix(100,100);
		for(int i = 0; i < 100; i++)
		{
			int x = (int) (Math.random() * 100);
			int y = (int) (Math.random() * 100);
			double v =  (Math.random() * 4000);
			largeDense.set(v,y,x);
			largeDense2.set(v,y,x);
		}
		assertTrue(largeDense.almostEqual(largeDense2));
		for(int i = 0; i < 1000; i++)
		{
			int x = (int) (Math.random() * 100);
			int y = (int) (Math.random() * 100);
			double v = (int) (Math.random() * 4000);
			largeDense.add(v,y,x);
			largeDense2.add(v,y,x);
			largeDense.add(v,y,x);
			largeDense2.add(v,y,x);
		}
		for(int i = 0; i < 100; i ++)
		{
			for(int j = 0; j < 100; j++)
				assertEquals(largeDense.at(i,j),largeDense2.at(i,j));
		}
		assertEquals(largeDense, largeDense2);
		SparseMatrix largeDense3 = (SparseMatrix)largeDense.add(largeDense);
		assertEquals(largeDense, largeDense2);
		SparseMatrix largeDense4 = largeDense.mul(0.1);
		SparseMatrix largeDenseRec =
			(SparseMatrix)largeDense.getLowerTriangleMatrix().add(largeDense.getDiagonalMatrix()).add(largeDense.getUpperTriangleMatrix());
		assertTrue(largeDense4.isSparse());
		assertEquals(largeDense, largeDense2);
		
		for(int i = 0; i < 100; i ++)
		{
			for(int j = 0; j < 100; j++)
				assertEquals(largeDense.at(i,j),largeDense2.at(i,j));
		}
		for(int i = 0; i < 100; i ++)
		{
			for(int j = 0; j < 100; j++)
			{
				assertTrue(Math.abs(largeDense4.at(i, j)- 0.1*largeDense.at(i, j)) < 1e-10);
				assertEquals(largeDenseRec.at(i, j), largeDense.at(i, j));
				assertEquals(largeDense3.at(i, j), 2*largeDense2.at(i, j));
			}
		}
		DenseVector largeVector = new DenseVector(100);
		SparseVector sparseVector = new SparseVector(100);
		for(int i = 0; i < 100; i++)
		{
			double x = Math.random();
			largeVector.set(x,i);
			sparseVector.set(x,i);
		}
		Vector mul1 = largeDense.mvMul(largeVector);
		Vector mul2 = largeDense3.mvMul(largeVector);
		Vector smul1 = largeDense.mvMul(sparseVector);
		Vector smul2 = largeDense3.mvMul(sparseVector);
		assertEquals(largeDense, largeDense3.mul(0.5));
		assertEquals(mul1.mul(2),mul2);
		assertTrue(mul1.almostEqual(mul2.mul(0.5),1e-14));
		for(int i = 0; i < 100; i++)
		{
			largeVector.add(Math.random(),i);
		}
		assertEquals(largeDense.transpose().transpose(), largeDense);
		assertTrue(largeDense.tvMul(largeVector).almostEqual(largeDense2.transpose().mvMul(largeVector),
			1e-12));
		Matrix m1 = largeDense.mmMul(largeDense2);
		Matrix m2 = largeDense2.mmMul(largeDense);
		assertTrue(m1.almostEqual(m2,1e-12));
	}
	@Test
	public void testSolve()
	{
		SparseMatrix largeDense = new SparseMatrix(50,50);
		DenseVector largeVector = new DenseVector(50);
		for(int i = 0; i < 4000; i++)
		{
			int x = (int) (Math.random() * 50);
			int y = (int) (Math.random() * 50);
			double v = (Math.random() * 4);
			double w = (Math.random() * 4);
			largeVector.set(w,x);
			largeDense.set(v,y,x);
			largeDense.add(5.2,i%50,i%50);
		}
	
		DenseMatrix denseInverse = largeDense.inverse();
		assertTrue(denseInverse.mvMul(largeVector).almostEqual(largeDense.solve(largeVector),1e-10));
		assertTrue(largeDense.mmMul(denseInverse).almostEqual(SparseMatrix.identity(50),1e-10));
		assertTrue(denseInverse.mmMul(largeDense).almostEqual(SparseMatrix.identity(50),1e-10));
		assertTrue(denseInverse.mmMul(largeDense).almostEqual(largeDense.mmMul(denseInverse),1e-10));
		assertTrue(denseInverse.inverse().almostEqual(largeDense,1e-10));
	
		assertTrue(new IterativeSolver<SparseMatrix>().solveBiCGStab(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		assertTrue(largeDense.mvMul(new IterativeSolver<SparseMatrix>().solveBiCGStab(largeDense,largeVector,
			1e-12)).almostEqual(largeVector,1e-8));
		assertTrue(new IterativeSolver<SparseMatrix>().solveGMRES(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		assertTrue(new IterativeSolver<SparseMatrix>().solveBiCGStab(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		assertTrue(new IterativeSolver<SparseMatrix>().solveGMRES(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
	
		SparseMatrix symmDense = (SparseMatrix)largeDense.add(largeDense.transpose());
		assertTrue(symmDense.solveSymm(largeVector).almostEqual(new IterativeSolver<SparseMatrix>().solveCG(symmDense,largeVector,1e-10),1e-8));
		assertTrue(symmDense.solve(largeVector).almostEqual(new IterativeSolver<SparseMatrix>().solveCG(symmDense,largeVector,1e-10),1e-8));
		assertTrue(symmDense.solve(largeVector).almostEqual(symmDense.solveSymm(largeVector),1e-10));
		
	}
	@Test
	public void testSpvecSolve()
	{
		SparseVector small = new SparseVector(40);
		small.add(4.3,4);
		small.add(7.8,45);
		assertEquals(small.euclidianNorm(), Math.sqrt(4.3*4.3+7.8*7.8));
		SparseMatrix largeDense = new SparseMatrix(500,500);
		SparseVector largeVector = new SparseVector(500);
		DenseVector largeVector2 = new DenseVector(500);
		for(int i = 0; i < 4000; i++)
		{
			int x = (int) (Math.random() * 500);
			int y = (int) (Math.random() * 500);
			double v = (Math.random() * 4);
			if(i%10 == 0)
			{
				double w = (Math.random() * 4);
				largeVector.set(w,x);
				largeVector2.set(w,x);
			}
			if(i%15 == 0)
			{
				double w = (Math.random() * 4);
				largeVector.add(w,x);
				largeVector2.add(w,x);
			}
			largeDense.set(v,y,x);
			largeDense.add(5.2,i%500,i%500);
		}

		assertTrue(largeVector.almostEqual(largeVector2));
		assertTrue(largeVector.add(largeVector).almostEqual(largeVector2.mul(2)));
		assertTrue(largeVector.sub(largeVector).euclidianNorm() <1e-14);
		assertTrue(largeVector2.sub(largeVector).euclidianNorm() <1e-14);
		assertTrue(largeVector2.sub(largeVector2).euclidianNorm() <1e-14);
		assertTrue(largeVector.sub(largeVector2).euclidianNorm() <1e-14);
		assertTrue(largeDense.mvMul(largeVector).almostEqual(largeDense.mvMul(largeVector)));
		assertTrue(largeDense.mvMul(largeVector).almostEqual(largeDense.mvMul(largeVector)));

		DenseVector sparsecopy = new DenseVector(largeVector);

		assertTrue(largeDense.mvMul(new IterativeSolver<SparseMatrix>().solveBiCGStab(largeDense,sparsecopy,
			1e-12)).almostEqual(sparsecopy,1e-8));

		assertTrue(new IterativeSolver<SparseMatrix>().solveBiCGStab(largeDense,largeVector,1e-12).almostEqual(new IterativeSolver<SparseMatrix>().solveGMRES(largeDense,largeVector,1e-10),
			1e-8));
		SparseMatrix symmDense = (SparseMatrix) largeDense.add(largeDense.transpose());

		assertTrue(symmDense.solveSymm(largeVector).almostEqual(new IterativeSolver<SparseMatrix>().solveCG(symmDense,largeVector,1e-10),1e-8));
		
	}
	
}
