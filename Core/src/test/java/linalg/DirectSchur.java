package linalg;

import scala.Function2;

class DirectSchur
	extends SchurSolver
{
	public DirectSchur(final BlockSparseMatrix blockMatrix)
	{
		super(blockMatrix);
	}
	
	@Override
	Function2<DenseMatrix, Vector, Vector> solveSchur()
	{
		return DenseMatrix::solve;
	}
}
