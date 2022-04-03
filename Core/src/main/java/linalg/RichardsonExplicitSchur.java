package linalg;

import basic.MetricWindow;
import io.vavr.Function2;

public class RichardsonExplicitSchur
	extends ExplicitSchurSolver
{
	public VectorMultiplyable preconditioner = null;
	
	public RichardsonExplicitSchur(final BlockSparseMatrix blockMatrix)
	{
		super(blockMatrix);
		System.out.println("inverting first" + blockMatrix.getBlockSizes()[0]);
	}
	
	@Override
	protected Function2<SparseMatrix, Vector, Vector> solveSchur()
	{
		return (A, b) ->
		{
			final IterativeSolverConvergenceMetric metric = new IterativeSolverConvergenceMetric(1e-7);
			MetricWindow.getInstance()
			            .setMetric("richSchur", metric);
			Vector iterate = b.mul(0);
			Vector residual = b.sub(A.mvMul(iterate));
			while (residual.euclidianNorm() > 1e-7)
			{
				iterate = iterate.add(preconditioner.mvMul(residual));
				residual = b.sub(A.mvMul(iterate));
				System.out.println();
				System.out.println();
				System.out.println("SCHURRES     " + residual.euclidianNorm());
				metric.publishIterate(residual.euclidianNorm());
				System.out.println();
				System.out.println();
			}
			return iterate;
		};
	}
}
