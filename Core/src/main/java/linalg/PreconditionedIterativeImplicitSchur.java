package linalg;

import scala.Function2;

public class PreconditionedIterativeImplicitSchur
	extends ImplicitSchurSolver
{
	final IterativeSolver it = new IterativeSolver();
	public VectorMultiplyable preconditioner;
	
	public PreconditionedIterativeImplicitSchur(final BlockSparseMatrix blockMatrix,
	                                            final VectorMultiplyable preconditioner)
	{
		super(blockMatrix);
		this.preconditioner = preconditioner;
		it.gm.MAX_RESTARTS = 0;
		it.showProgress = true;
	}
	
	@Override
	protected Function2<VectorMultiplyable, Vector, Vector> solveSchur()
	{
		return (A, b) ->
		{
			it.gm.ITERATIONS_BEFORE_RESTART = 40;
			System.out.println(A.getVectorSize());
//			System.out.println(((AMGPreconditionerSpace) preconditioner).finest_system.mvMul(b)
//			                                                                          .sub(A.mvMul(b))
//			                                                                          .absMaxElement());
			return it.solvePGMRES(A, preconditioner, b, 1e-7);
		};
	}
}
