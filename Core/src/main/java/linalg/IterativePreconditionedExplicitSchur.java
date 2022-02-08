package linalg;

import scala.Function2;

public class IterativePreconditionedExplicitSchur
	extends ExplicitSchurSolver
{
	final IterativeSolver it = new IterativeSolver();
	public VectorMultiplyable preconditioner = null;
	
	public IterativePreconditionedExplicitSchur(final BlockSparseMatrix blockMatrix)
	{
		super(blockMatrix);
		System.out.println("inverting first" + blockMatrix.getBlockSizes()[0]);
	}
	
	@Override
	protected Function2<SparseMatrix, Vector, Vector> solveSchur()
	{
		return (A, b) ->
		{
			if (preconditioner == null)
				return it.solveGMRES(A, b, 1e-9);
			else
				return it.solvePGMRES(A, preconditioner, b, 1e-6);
		};
	}
}
