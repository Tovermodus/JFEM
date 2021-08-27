package linalg;

public interface DirectlySolvable extends Matrix
{
	default DirectlySolvable inverse()
	{
		return (new DenseMatrix(this)).inverse();
	}
	
	Vector solve(Vector rhs);
	
	Vector solveSymm(Vector rhs);
}
