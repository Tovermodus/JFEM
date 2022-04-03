package linalg;

import io.vavr.Function2;

public class IterativeExplicitSchur
	extends ExplicitSchurSolver
{
	final IterativeSolver it = new IterativeSolver();
	DenseMatrix firstInverse;
	
	public IterativeExplicitSchur(final BlockSparseMatrix blockMatrix)
	{
		super(blockMatrix);
		System.out.println("inverting first" + blockMatrix.getBlockSizes()[0]);
		firstInverse = blockMatrix.getBlockMatrix(0, 0)
		                          .inverse();
	}
	
	@Override
	protected Function2<SparseMatrix, Vector, Vector> solveSchur()
	{
		return (A, b) -> it.solvePGMRES(A, firstInverse, b, 1e-9);
	}
}
