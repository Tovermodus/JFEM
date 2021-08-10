package linalg;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class PowerIterationTest
{
	@Test
	public void testSmall()
	{
		DenseMatrix d = DenseMatrix.squareMatrixFromValues(2,1,0,1,2,0,0,0,1);
		assertTrue(Math.abs(d.powerIterationSymmetric() - 3) < 1e-5);
	}
	@Test
	public void testIdentity()
	{
		SparseMatrix s = SparseMatrix.identity(1000);
		assertTrue(Math.abs(s.powerIterationSymmetric() - 1) < 1e-5);
		
	}
	@Test
	public void testToeplitz()
	{
		int n = 100;
		SparseMatrix s = SparseMatrix.identity(n).mul(4);
		for(int i = 1 ; i < n-1; i++)
		{
			s.add(5,i,i+1);
			s.add(5,i,i-1);
		}
		s.add(5,0,1);
		s.add(5,n-1,n-2);
		double eigenVal = 4 + 2 * 5 * Math.cos(1 * Math.PI / (n -1));
		assertTrue(Math.abs(s.powerIterationSymmetric() - eigenVal) < 1e-5);
	}
	
}
