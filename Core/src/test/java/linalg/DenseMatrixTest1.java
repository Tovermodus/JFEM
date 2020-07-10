package linalg;

import com.google.common.primitives.Doubles;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;
public class DenseMatrixTest1
{

	@Test
	public void testGetters() {
		DenseVector vector;
		DenseMatrix singularMatrix;
		DenseMatrix invertibleMatrix;
		vector = DenseVector.vectorFromValues(0,1,2,3);
		singularMatrix = DenseMatrix.squareMatrixFromValues(0,1,0,5);
		invertibleMatrix = DenseMatrix.squareMatrixFromValues(0,1,1,1);
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
		assertEquals(vector.euclidianNorm(), Math.sqrt(1+4+9));
		
	}
	@Test
	public void testLarge()
	{
		DenseMatrix largeDense = new DenseMatrix(50,50);
		DenseMatrix largeDense2 = new DenseMatrix(50,50);
		for(int i = 0; i < 100; i++)
		{
			int x = (int) (Math.random() * 50);
			int y = (int) (Math.random() * 50);
			double v =  (Math.random() * 4000);
			largeDense.set(v,y,x);
			largeDense2.set(v,y,x);
		}
		for(int i = 0; i < 50; i ++)
		{
			for(int j = 0; j < 50; j++)
				assertEquals(largeDense.at(i,j),largeDense2.at(i,j));
		}
		assertTrue(largeDense.almostEqual(largeDense2));
		for(int i = 0; i < 1000; i++)
		{
			int x = (int) (Math.random() * 50);
			int y = (int) (Math.random() * 50);
			double v = (int) (Math.random() * 4000);
			largeDense.add(v,y,x);
			largeDense2.add(v,y,x);
		}
		assertEquals(largeDense, largeDense2);
		DenseMatrix largeDense3 = largeDense.add(largeDense);
		assertEquals(largeDense, largeDense2);
		DenseMatrix largeDense4 = largeDense.mul(0.1);
		DenseMatrix largeDenseRec =
			largeDense.getLowerTriangleMatrix().add(largeDense.getDiagonalMatrix()).add(largeDense.getUpperTriangleMatrix());
		assertFalse(largeDense4.isSparse());
		assertEquals(largeDense, largeDense2);
		for(int i = 0; i < 50; i ++)
		{
			for(int j = 0; j < 50; j++)
			{
				assertEquals(largeDense3.at(i, j), 2*largeDense.at(i, j));
				assertEquals(largeDense4.at(i, j), 0.1*largeDense.at(i, j));
				assertEquals(largeDenseRec.at(i, j), largeDense.at(i, j));
			}
		}
		DenseVector largeVector = new DenseVector(50);
		SparseVector sparseVector = new SparseVector(50);
		for(int i = 0; i < 50; i++)
		{
			double x = Math.random();
			largeVector.set(x,i);
			sparseVector.set(x,i);
		}
		DenseVector mul1 = largeDense.mvMul(largeVector);
		DenseVector mul2 = largeDense3.mvMul(largeVector);
		DenseVector smul1 = largeDense.mvMul(sparseVector);
		DenseVector smul2 = largeDense3.mvMul(sparseVector);
		assertEquals(largeDense, largeDense3.mul(0.5));
		for(int i = 0; i < 50; i++)
		{
			assertEquals(mul1.mul(2).at(i), mul2.at(i));
			assertEquals(mul1.at(i)*2, mul2.at(i));
		}
		assertTrue(mul1.almostEqual(mul2.mul(0.5),1e-14));
		for(int i = 0; i < 50; i++)
		{
			largeVector.add(Math.random(),i);
		}
		assertEquals(largeDense.transpose().transpose(), largeDense);
		assertTrue(largeDense.tvMul(largeVector).almostEqual(largeDense2.transpose().mvMul(largeVector),
			1e-12));
	}
	@Test
	public void testSolve()
	{
		DenseMatrix largeDense = new DenseMatrix(50,50);
		DenseVector largeVector = new DenseVector(50);
		for(int i = 0; i < 100; i++)
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
		assertTrue(denseInverse.mmMul(largeDense).almostEqual(largeDense.mmMul(denseInverse),1e-10));
		assertTrue(denseInverse.mmMul(largeDense).almostEqual(DenseMatrix.identity(50),1e-10));
		assertTrue(denseInverse.inverse().almostEqual(largeDense,1e-10));
		assertTrue(new IterativeSolver<DenseMatrix>().solveBiCGStab(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		assertTrue(largeDense.mvMul(new IterativeSolver<DenseMatrix>().solveBiCGStab(largeDense,largeVector,
			1e-12)).almostEqual(largeVector,1e-8));
		assertTrue(new IterativeSolver<DenseMatrix>().solveGMRES(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		assertTrue(new IterativeSolver<DenseMatrix>().solveBiCGStab(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		assertTrue(new IterativeSolver<DenseMatrix>().solveGMRES(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		DenseMatrix symmDense = denseInverse.add(denseInverse.transpose());
		assertTrue(symmDense.solveSymm(largeVector).almostEqual(new IterativeSolver<DenseMatrix>().solveCG(symmDense,largeVector,1e-10),1e-8));
		assertTrue(symmDense.solve(largeVector).almostEqual(new IterativeSolver<DenseMatrix>().solveCG(symmDense,largeVector,1e-10),1e-8));
		assertTrue(symmDense.solve(largeVector).almostEqual(symmDense.solveSymm(largeVector),1e-10));
		
	}
	@Test
	public void testSpvecSolve()
	{
		SparseVector small = new SparseVector(40);
		small.add(4.3,4);
		small.add(7.8,45);
		assertEquals(small.euclidianNorm(), Math.sqrt(4.3*4.3+7.8*7.8));
		DenseMatrix largeDense = new DenseMatrix(50,50);
		SparseVector largeVector = new SparseVector(50);
		DenseVector largeVector2 = new DenseVector(50);
		for(int i = 0; i < 100; i++)
		{
			int x = (int) (Math.random() * 50);
			int y = (int) (Math.random() * 50);
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
			largeDense.add(5.2,i%50,i%50);
		}
		DenseMatrix denseInverse = largeDense.inverse();
		for(int j = 0; j < 50; j++)
		{
			assertEquals(largeVector.at(j), largeVector2.at(j));
		}
		for(int j = 0; j < 50; j++)
		{
			assertEquals(largeVector.mul(0.1).at(j), largeVector2.at(j)*0.1);
		}
		assertTrue(largeVector.almostEqual(largeVector2));
		assertTrue(largeVector.add(largeVector).almostEqual(largeVector2.mul(2)));
		assertTrue(largeVector.sub(largeVector).euclidianNorm() <1e-14);
		assertTrue(largeVector2.sub(largeVector).euclidianNorm() <1e-14);
		assertTrue(largeVector2.sub(largeVector2).euclidianNorm() <1e-14);
		assertTrue(largeVector.sub(largeVector2).euclidianNorm() <1e-14);
		assertTrue(largeDense.mvMul(largeVector).almostEqual(largeDense.mvMul(largeVector)));
		assertTrue(largeDense.mvMul(largeVector).almostEqual(largeDense.mvMul(largeVector)));
		assertTrue(denseInverse.mvMul(largeVector).almostEqual(largeDense.solve(largeVector),1e-10));
		assertTrue(denseInverse.mmMul(largeDense).almostEqual(largeDense.mmMul(denseInverse),1e-10));
		assertTrue(denseInverse.mmMul(largeDense).almostEqual(DenseMatrix.identity(50),1e-10));
		assertTrue(denseInverse.inverse().almostEqual(largeDense,1e-10));
		denseInverse = largeDense.inverse();
		DenseVector sparsecopy = new DenseVector(largeVector);
		assertTrue(largeDense.mvMul(new IterativeSolver<DenseMatrix>().solveBiCGStab(largeDense,sparsecopy,
			1e-12)).almostEqual(sparsecopy,1e-8));
		assertTrue(largeDense.mvMul(new IterativeSolver<DenseMatrix>().solveBiCGStab(largeDense,largeVector,
			1e-12)).almostEqual(sparsecopy,1e-8));
		assertTrue(new IterativeSolver<DenseMatrix>().solveBiCGStab(largeDense,largeVector,1e-12).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		assertTrue(new IterativeSolver<DenseMatrix>().solveGMRES(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		assertTrue(new IterativeSolver<DenseMatrix>().solveBiCGStab(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		assertTrue(new IterativeSolver<DenseMatrix>().solveGMRES(largeDense,largeVector,1e-10).almostEqual(denseInverse.mvMul(largeVector),
			1e-8));
		DenseMatrix symmDense = denseInverse.add(denseInverse.transpose());
		assertTrue(symmDense.solveSymm(largeVector).almostEqual(new IterativeSolver<DenseMatrix>().solveCG(symmDense,largeVector,1e-10),1e-8));
		assertTrue(symmDense.solve(largeVector).almostEqual(new IterativeSolver<DenseMatrix>().solveCG(symmDense,largeVector,1e-10),1e-8));
		assertTrue(symmDense.solve(largeVector).almostEqual(symmDense.solveSymm(largeVector),1e-10));
		
	}
	
}
