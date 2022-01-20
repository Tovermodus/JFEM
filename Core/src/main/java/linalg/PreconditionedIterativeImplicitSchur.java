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
		this.preconditioner = //preconditioner;
			new BlockDenseMatrix(blockMatrix.getBlockMatrix(0, 0),
			                     blockMatrix.getBlockSizes()[0] / 100)
				.getInvertedDiagonalMatrix();
		it.gm.MAX_RESTARTS = 0;
		System.out.println("inverting first" + blockMatrix.getBlockSizes()[0]);
	}
	
	@Override
	Function2<VectorMultiplyable, Vector, Vector> solveSchur()
	{
		return (A, b) ->
		{
//			Vector iterate = new DenseVector(b);
//			for (int i = 0; i < 40; i++)
//				iterate = preconditioner.mvMul(iterate);
//			return iterate;
			return it.solvePGMRES(A, preconditioner, b, 1e-10);
		};
	}
}
