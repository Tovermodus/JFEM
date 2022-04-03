package linalg;

import io.vavr.Function2;

public class IterativePreconditionedExplicitSchur
	extends ExplicitSchurSolver
{
	GMRES2 gm = new GMRES2(1e-6, true);
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
				return it.solvePCG(A, preconditioner, b, 1e-6);
		};
	}
}
