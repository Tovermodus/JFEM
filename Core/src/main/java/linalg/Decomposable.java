package linalg;

public interface Decomposable
	extends Matrix
{
	Matrix getStrictlyLowerTriangleMatrix();
	
	Matrix getStrictlyUpperTriangleMatrix();
	
	Matrix getDiagonalMatrix();
}
