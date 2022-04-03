package dlm;

import basic.MetricWindow;
import io.vavr.Function2;
import linalg.*;
import multigrid.MGInterface;

public class MGSchur
	extends ImplicitSchurSolver
{
	final IterativeSolver it = new IterativeSolver();
	final MGInterface preconditioner;
	DenseVector iterate;
	
	public MGSchur(final BlockSparseMatrix blockMatrix,
	               final MGInterface preconditioner)
	{
		super(blockMatrix);
		this.preconditioner = preconditioner;
//			new BlockDenseMatrix(blockMatrix.getBlockMatrix(0, 0),
//			                     blockMatrix.getBlockSizes()[0] / 100)
//				.getInvertedDiagonalMatrix();
		it.gm.MAX_RESTARTS = 0;
		System.out.println("inverting first" + blockMatrix.getBlockSizes()[0]);
	}
	
	@Override
	protected Function2<VectorMultiplyable, Vector, Vector> solveSchur()
	{
		return (A, b) ->
		{
			final IterativeSolverConvergenceMetric icm = new IterativeSolverConvergenceMetric(1e-8);
			MetricWindow.getInstance()
			            .setMetric("schurmg", icm);
			iterate = new DenseVector(preconditioner.getFinestSystem()
			                                        .getVectorSize());
			for (int i = 0; i < 100; i++)
			{
				iterate = new DenseVector(preconditioner.vCycle(iterate, b));
				final double res = preconditioner.getFinestSystem()
				                                 .mvMul(iterate)
				                                 .sub(b)
				                                 .absMaxElement();
				System.out.println("mgcycle " + i + " " + res);
				icm.publishIterate(res);
//				System.out.println(iterate.absMaxElement());
			}
			return b;//iterate;
		};
	}
}
