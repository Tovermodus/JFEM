package linalg;

import scala.Function2;

public class IterativeImplicitSchur
	extends ImplicitSchurSolver
{
	final IterativeSolver it = new IterativeSolver();
	DenseMatrix firstInverse;
	
	public IterativeImplicitSchur(final BlockSparseMatrix blockMatrix)
	{
		super(blockMatrix);
		System.out.println("inverting first" + blockMatrix.getBlockSizes()[0]);
		firstInverse = blockMatrix.getBlockMatrix(0, 0)
		                          .inverse();
	}
	
	@Override
	protected Function2<VectorMultiplyable, Vector, Vector> solveSchur()
	{
		return (A, b) -> it.solvePGMRES(A, firstInverse, b, 1e-9);
	}
}
