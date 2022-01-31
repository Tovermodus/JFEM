package linalg;

import scala.Function2;

public class PreconditionedIterativeImplicitSchur
	extends ImplicitSchurSolver
{
	final IterativeSolver it = new IterativeSolver();
	final VectorMultiplyable preconditioner;
	
	public PreconditionedIterativeImplicitSchur(final BlockSparseMatrix blockMatrix,
	                                            final VectorMultiplyable preconditioner)
	{
		super(blockMatrix);
		this.preconditioner = preconditioner;
		it.gm.MAX_RESTARTS = 0;
		it.showProgress = true;
		System.out.println("inverting first" + blockMatrix.getBlockSizes()[0]);
	}
	
	@Override
	protected Function2<VectorMultiplyable, Vector, Vector> solveSchur()
	{
		return (A, b) ->
		{
			it.gm.ITERATIONS_BEFORE_RESTART = 50;
			System.out.println(A.getVectorSize());
			return it.solvePGMRES(A, preconditioner, b, 1e-10);
			//new GMRES2(1e-10, true).solve(A, preconditioner, b);
		};
	}
}
