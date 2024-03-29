package linalg;

import basic.MetricWindow;
import basic.StatLogger;
import io.vavr.Function2;
import multigrid.MGPreconditionerInterface;

public class MultigridFixedpointSchur
	extends ImplicitSchurSolver
{
	public MGPreconditionerInterface<?, ?, ?, ?, ?, ?, ?> mg;
	
	public MultigridFixedpointSchur(final BlockSparseMatrix blockMatrix,
	                                final MGPreconditionerInterface<?, ?, ?, ?, ?, ?, ?> preconditioner)
	{
		super(blockMatrix);
		this.mg = preconditioner;
	}
	
	@Override
	protected Function2<VectorMultiplyable, Vector, Vector> solveSchur()
	{
		return (A, b) ->
		{
			DenseVector iterate = new DenseVector(b);
			final double initialres = b.sub(A.mvMul(iterate))
			                           .euclidianNorm();
			double res = initialres;
			int it;
			final IterativeSolverConvergenceMetric cm = new IterativeSolverConvergenceMetric(1e-8);
			MetricWindow.getInstance()
			            .setMetric("mg", cm);
			for (it = 0; it < 100 && res > 1e-6; it++)
			
			{
				cm.publishIterate(res);
				iterate = new DenseVector(mg.vCycle(iterate, b));
				res = mg.getFinestSystem()
				        .mvMul(iterate)
				        .sub(b)
				        .euclidianNorm();
				System.out.println(res / initialres + " after iterations " + it);
				StatLogger.log("MGFIxedpointSchur: " + res + " iteration " + it + " initialres " + initialres);
			}
			return iterate;
		};
	}
}
