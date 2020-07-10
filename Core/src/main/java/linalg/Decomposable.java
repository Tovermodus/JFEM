package linalg;

import linalg.Matrix;

public interface Decomposable extends Matrix
{
	Matrix getLowerTriangleMatrix();
	Matrix getUpperTriangleMatrix();
	Matrix getDiagonalMatrix();
}
