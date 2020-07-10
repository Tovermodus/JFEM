package linalg;

import linalg.Matrix;
import linalg.Vector;

public interface DirectlySolvable extends Matrix
{
	Matrix inverse();
	Vector solve(Vector rhs);
	Vector solveSymm(Vector rhs);
}
