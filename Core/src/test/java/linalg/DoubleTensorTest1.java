package linalg;

import com.google.common.primitives.Doubles;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;
public class DoubleTensorTest1
{
	
	@Test
	public void testGetters() {
		DoubleTensor vector;
		DoubleTensor singularMatrix;
		DoubleTensor invertibleMatrix;
		DoubleTensor singularSparseMatrix;
		DoubleTensor invertibleSparseMatrix;
		vector = DoubleTensor.vectorFromValues(0,1,2,3);
		singularMatrix = DoubleTensor.squareMatrixFromValues(0,1,0,5);
		invertibleMatrix = DoubleTensor.squareMatrixFromValues(0,1,1,1);
		singularSparseMatrix = new DoubleTensor(3,4,true);
		singularSparseMatrix.add(0,0,1);
		singularSparseMatrix.add(1,1,-1);
		singularSparseMatrix.add(1,2,2);
		singularSparseMatrix.add(2,1,2);
		singularSparseMatrix.add(2,2,-4);
		invertibleSparseMatrix = new DoubleTensor(3,3,true);
		invertibleSparseMatrix.add(0,0,3);
		invertibleSparseMatrix.add(1,1,2);
		invertibleSparseMatrix.add(0,2,-1);
		invertibleSparseMatrix.add(2,2,99);
		assertEquals(1, vector.getDimension());
		assertEquals(2, singularMatrix.getDimension());
		assertEquals(2, invertibleMatrix.getDimension());
		assertEquals(2, singularSparseMatrix.getDimension());
		assertEquals(2, invertibleSparseMatrix.getDimension());
		assertEquals(1, vector.getCols());
		assertEquals(2, singularMatrix.getCols());
		assertEquals(2, invertibleMatrix.getCols());
		assertEquals(4, singularSparseMatrix.getCols());
		assertEquals(3, invertibleSparseMatrix.getCols());
		assertEquals(4, vector.getRows());
		assertEquals(2, singularMatrix.getRows());
		assertEquals(2, invertibleMatrix.getRows());
		assertEquals(3, singularSparseMatrix.getRows());
		assertEquals(3, invertibleSparseMatrix.getRows());
		assertFalse(vector.isSparse());
		assertFalse(singularMatrix.isSparse());
		assertFalse(invertibleMatrix.isSparse());
		assertTrue(singularSparseMatrix.isSparse());
		assertTrue(invertibleSparseMatrix.isSparse());
		assertEquals(4, vector.size());
		assertEquals(4, singularMatrix.size());
		assertEquals(4, invertibleMatrix.size());
		assertEquals(12, singularSparseMatrix.size());
		assertEquals(9, invertibleSparseMatrix.size());
		assertEquals(vector.vectorNorm(), Math.sqrt(1+4+9));
	}
	@Test
	public void testLarge()
	{
		DoubleTensor largeDense = new DoubleTensor(50,50,false);
		DoubleTensor largeSparse = new DoubleTensor(50,50,true);
		for(int i = 0; i < 100; i++)
		{
			int x = (int) (Math.random() * 50);
			int y = (int) (Math.random() * 50);
			int v = (int) (Math.random() * 4000);
			largeDense.set(y,x,v);
			largeSparse.set(y,x,v);
		}
		for(int i = 0; i < 50; i ++)
		{
			for(int j = 0; j < 50; j++)
				assertEquals(largeDense.at(i,j),largeSparse.at(i,j));
		}
		for(int i = 0; i < 1000; i++)
		{
			int x = (int) (Math.random() * 50);
			int y = (int) (Math.random() * 50);
			int v = (int) (Math.random() * 4000);
			largeDense.add(y,x,v);
			largeSparse.add(y,x,v);
		}
		assertEquals(largeDense, largeSparse);
		DoubleTensor largeSparse2 = largeSparse.add(largeSparse);
		assertEquals(largeDense, largeSparse);
		DoubleTensor largeDense2 = largeDense.add(largeSparse);
		assertEquals(largeDense, largeSparse);
		DoubleTensor largeDense3 = largeDense.add(largeDense);
		DoubleTensor largeDense4 = largeDense.mul(0.1);
		DoubleTensor largeDenseRec =
			largeDense.getLowerTriangleMatrix().add(largeDense.getDiagonalMatrix()).add(largeDense.getUpperTriangleMatrix());
		assertTrue(largeSparse2.isSparse());
		assertFalse(largeDense2.isSparse());
		assertEquals(largeDense, largeSparse);
		for(int i = 0; i < 50; i ++)
		{
			for(int j = 0; j < 50; j++)
			{
				assertEquals(2*largeSparse.at(i, j), 2*largeDense.at(i, j));
				assertEquals(largeSparse2.at(i, j), 2*largeSparse.at(i, j));
				assertEquals(largeDense2.at(i, j), 2*largeDense.at(i, j));
				assertEquals(largeDense2.at(i, j), largeDense3.at(i, j));
				assertEquals(largeDense2.at(i, j), largeSparse2.at(i, j));
				assertEquals(largeDense4.at(i, j), 0.1*largeSparse.at(i, j));
				assertEquals(largeDenseRec.at(i, j), largeSparse.at(i, j));
			}
		}
		DoubleTensor largeVector = new DoubleTensor(50);
		for(int i = 0; i < 50; i++)
		{
			largeVector.add(i,Math.random());
		}
		DoubleTensor mul1 = largeDense.mvmul(largeVector);
		DoubleTensor mul2 = largeDense2.mvmul(largeVector);
		DoubleTensor mul3 = largeSparse.mvmul(largeVector);
		assertEquals(largeDense, largeSparse);
		for(int i = 0; i < 50; i++)
		{
			assertEquals(mul1.mul(2).at(i), mul2.at(i));
		}
		assertTrue(mul1.almostequals(mul3,1e-14));
		for(int i = 0; i < 50; i++)
		{
			largeVector.add(i,Math.random());
		}
		mul1 = largeDense.mvmul(largeVector);
		mul2 = largeDense2.mvmul(largeVector);
		mul3 = largeSparse.mvmul(largeVector);
		assertTrue(mul1.almostequals(mul3,1e-14));
		for(int i = 0; i < 50; i++)
		{
			assertEquals(mul1.mul(2).at(i), mul2.at(i));
		}
		assertEquals(largeDense.transpose().transpose(), largeDense);
		assertEquals(largeSparse.transpose().transpose(), largeSparse);
		assertTrue(largeDense.tvmul(largeVector).almostequals(largeSparse.transpose().mvmul(largeVector),
			1e-12));
		assertTrue(largeSparse.tvmul(largeVector).almostequals(largeSparse.transpose().mvmul(largeVector),
			1e-12));
	}
	@Test
	public void testSolve()
	{
		DoubleTensor largeDense = new DoubleTensor(50,50,false);
		DoubleTensor largeVector = new DoubleTensor(50);
		DoubleTensor largeSparse = new DoubleTensor(50,50,true);
		for(int i = 0; i < 100; i++)
		{
			int x = (int) (Math.random() * 50);
			int y = (int) (Math.random() * 50);
			double v = (Math.random() * 4);
			double w = (Math.random() * 4);
			largeVector.set(x,w);
			largeDense.set(y,x,v);
			largeSparse.set(y,x,v);
			largeDense.add(i%50,i%50,5.2);
			largeSparse.add(i%50,i%50,5.2);
		}
		DoubleTensor denseInverse = largeDense.inverse();
		DoubleTensor sparseInverse = largeDense.inverse();
		assertTrue(denseInverse.almostequals(sparseInverse,1e-8));
		assertTrue(denseInverse.mvmul(largeVector).almostequals(largeDense.solve(largeVector),1e-10));
		assertTrue(sparseInverse.mvmul(largeVector).almostequals(largeSparse.solve(largeVector),1e-10));
		assertTrue(sparseInverse.mvmul(largeVector).almostequals(largeDense.solve(largeVector),1e-8));
		assertTrue(sparseInverse.mmmul(largeDense).almostequals(denseInverse.mmmul(largeDense),1e-10));
		assertTrue(sparseInverse.mmmul(largeDense).almostequals(largeSparse.mmmul(sparseInverse),1e-10));
		assertTrue(sparseInverse.mmmul(largeSparse).almostequals(DoubleTensor.identity(50),1e-10));
		assertTrue(sparseInverse.inverse().almostequals(largeSparse,1e-10));
		assertTrue(sparseInverse.inverse().mvmul(largeVector).almostequals(largeSparse.mvmul(largeVector),
			1e-10));
		assertTrue(largeDense.solveBiCGStab(largeVector,1e-10).almostequals(denseInverse.mvmul(largeVector),
			1e-8));
		assertTrue(largeDense.solveGMRES(largeVector,1e-10).almostequals(denseInverse.mvmul(largeVector),
			1e-8));
		assertTrue(largeSparse.solveBiCGStab(largeVector,1e-10).almostequals(denseInverse.mvmul(largeVector),
			1e-8));
		assertTrue(largeSparse.solveGMRES(largeVector,1e-10).almostequals(denseInverse.mvmul(largeVector),
			1e-8));
		DoubleTensor symmDense = denseInverse.add(denseInverse.transpose());
		DoubleTensor symmSparse = denseInverse.add(denseInverse.transpose());
		assertEquals(symmDense, symmSparse);
		assertTrue(symmDense.solveSymm(largeVector).almostequals(symmSparse.solveCG(largeVector,1e-10),1e-8));
		assertTrue(symmDense.solve(largeVector).almostequals(symmDense.solveCG(largeVector,1e-10),1e-8));
		assertTrue(symmDense.solve(largeVector).almostequals(symmDense.solveSymm(largeVector),1e-10));
		
	}
	
}
