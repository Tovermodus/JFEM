package linalg;

import basic.MetricWindow;

public class GMRES2
{
	IterativeSolverConvergenceMetric icm;
	final double tol;
	public boolean verbose;
	Vector iterate;
	
	public GMRES2(final double tol)
	{
		this(tol, false);
	}
	
	public GMRES2(final double tol, final boolean show)
	{
		this.tol = tol;
		icm = new IterativeSolverConvergenceMetric(tol);
		if (show)
			MetricWindow.getInstance()
			            .setMetric("gmres2" + tol, icm);
	}
	
	public GMRES2(final String name, final double tol)
	{
		this.tol = tol;
		icm = new IterativeSolverConvergenceMetric(tol);
		MetricWindow.getInstance()
		            .setMetric(name, icm);
	}
	
	public Vector solve(final VectorMultiplyable A,
	                    final Vector b)
	{
		if (iterate == null)
		{
			iterate = new DenseVector(b.getLength());
			((DenseVector) iterate).set(1, 0);
		}
		return solveInternal(A, SparseMatrix.identity(b.getLength()), b, iterate);
	}
	
	public Vector solve(final VectorMultiplyable A,
	                    final VectorMultiplyable preconditioner,
	                    final Vector b)
	{
		if (iterate == null)
		{
			iterate = new DenseVector(b.getLength());
			((DenseVector) iterate).set(1, 0);
		}
		return solveInternal(A,
		                     preconditioner,
		                     b,
		                     iterate);
	}
	
	public Vector solveInternal(final VectorMultiplyable A,
	                            final VectorMultiplyable M,
	                            final Vector b,
	                            final Vector initial)
	{
		iterate = initial;
		Vector residual = b.sub(A.mvMul(iterate));
		Vector Mres = M.mvMul(residual);
		final VectorMultiplyable MA = VectorMultiplyable.concatenate(M, A);
		final Arnoldi a = new Arnoldi(Mres, MA);
		while (residual.euclidianNorm() > tol)
		{
			final Vector v = a.latestV();
			final Vector Av = a.latestAv();
			final double t = Mres.inner(Av) / Av.inner(Av);
			final Vector step = v.mul(t);
			iterate = iterate.add(step);
//			final Vector res2 = residual.sub(A.mvMul(v)
//			                                  .mul(t));
			residual = residual.sub(A.mvMul(v)
			                         .mul(t));// b.sub(A.mvMul(iterate));
//			System.out.println("GMRES TEESSST " + residual.sub(res2)
//			                                              .euclidianNorm());
			//final Vector Mres2 = Mres.sub(Av.mul(t));
			Mres = Mres.sub(Av.mul(t));//M.mvMul(residual);
//			System.out.println("GMRES TEESSST2 " + Mres.sub(Mres2)
//			                                           .euclidianNorm());
//			Vector MAV2 = M.mvMul(A.mvMul(v));
//			MAV2 = MAV2.mul(1. / MAV2.euclidianNorm());
//			System.out.println("GMRES TEESSST3 " + Av.sub(MAV2)
//			                                         .euclidianNorm());
			if (verbose)
			{
				System.out.println("residuaaaaaaal " + residual.euclidianNorm() + " MRES " + Mres.euclidianNorm());
			}
			icm.publishIterate(residual.euclidianNorm());
			if (a.arnoldiStep())
			{
				if (verbose)
					System.out.println("restart");
				return solveInternal(A, M, b, iterate);
			}
		}
		return iterate;
	}
}
