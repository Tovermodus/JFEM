package linalg;

import linalg.Matrix;
import linalg.Vector;

public interface DirectlySolvable extends Matrix
{
	default DenseMatrix inverse()
	{
		return (new DenseMatrix(this)).inverse();
	}
	Vector solve(Vector rhs);
	Vector solveSymm(Vector rhs);
}
