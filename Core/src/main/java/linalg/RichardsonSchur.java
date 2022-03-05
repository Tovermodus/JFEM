package linalg;

import basic.MetricWindow;
import scala.Function2;

public class RichardsonSchur
	extends ImplicitSchurSolver
{
	public VectorMultiplyable preconditioner;
	
	public RichardsonSchur(final BlockSparseMatrix blockMatrix,
	                       final VectorMultiplyable preconditioner)
	{
		super(blockMatrix);
		this.preconditioner = preconditioner;
	}
	
	@Override
	protected Function2<VectorMultiplyable, Vector, Vector> solveSchur()
	{
		return (A, b) ->
		{
			final IterativeSolverConvergenceMetric it = new IterativeSolverConvergenceMetric(1e-7);
			MetricWindow.getInstance()
			            .setMetric("schurrich", it);
			Vector iterate = b.mul(0);
			Vector residual = b.sub(A.mvMul(iterate));
			while (residual.euclidianNorm() > 1e-7)
			{
				iterate = iterate.add(preconditioner.mvMul(residual));
				residual = b.sub(A.mvMul(iterate));
				System.out.println();
				System.out.println();
				it.publishIterate(residual.euclidianNorm());
				System.out.println("SCHURRES     " + residual.euclidianNorm());
				System.out.println();
				System.out.println();
			}
			return iterate;
		};
	}
}
