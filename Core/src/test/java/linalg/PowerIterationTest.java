package linalg;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class PowerIterationTest
{
	@Test
	public void testSmall()
	{
		DenseMatrix d = DenseMatrix.squareMatrixFromValues(2,1,0,1,2,0,0,0,1);
		assertTrue(Math.abs(d.powerIteration() - 3) < 1e-5);
	}
	@Test
	public void testIdentity()
	{
		SparseMatrix s = SparseMatrix.identity(1000);
		assertTrue(Math.abs(s.powerIteration() - 1) < 1e-5);
		
	}
	@Test
	public void testToeplitz()
	{
		SparseMatrix s = SparseMatrix.identity(1000).mul(4);
		s.set(0,1,7);
		s.set(1000,999,5);
		for(int i = 1 ; i < 999; i++)
		{
			s.set(i,i+1,7);
			s.set(i,i-1,5);
		}
		assertTrue(Math.abs(s.powerIteration() - 4) < 1e-5);
	}
	
}
