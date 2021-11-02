package linalg;

import basic.Interruptor;
import basic.MetricWindow;
import basic.PerformanceArguments;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class IterativeSolver
{
	private GMRES gm;
	ExecutorService ex;
	public boolean showProgress = false;
	public boolean showInterrupt = false;
	Vector lastSolution = null;
	public IterativeSolverConvergenceMetric metric;
	
	public IterativeSolver()
	{
		metric = new IterativeSolverConvergenceMetric(1);
		MetricWindow.getInstance()
		            .addMetric(metric);
	}
	
	public IterativeSolver(final boolean showProgress)
	{
		metric = new IterativeSolverConvergenceMetric(1);
		MetricWindow.getInstance()
		            .addMetric(metric);
	}
	
	public Vector solveCG(final VectorMultiplyable operator, final Vector rhs, final double tol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			final DenseVector d = new DenseVector(rhs.getLength());
			for (int i = 0; i < d.getLength(); i++)
				d.set(Math.random(), i);
			if (operator.mvMul(d)
			            .absMaxElement() < 1e-14)
				System.err.println("operator has a kernel");
		}
		metric.goal = tol;
		metric.restart();
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showInterrupt)
			ex.execute(i);
		final int n = rhs.getLength();
		Vector z;
		Vector newResiduum;
		DenseVector iterate = new DenseVector(n);
		iterate.set(1, 0);
		if (lastSolution != null)
			iterate = new DenseVector(lastSolution);
		final Vector m = operator.mvMul(iterate);
		Vector residuum = rhs.sub(m);
		Vector defect = new DenseVector(residuum);
		double alpha;
		double beta;
		for (int iter = 0; iter < n && residuum.euclidianNorm() > tol && i.running; iter++)
		{
			z = operator.mvMul(defect);
			alpha = residuum.inner(residuum) / defect.inner(z);
			if (alpha < 0 && showProgress)
				System.out.println("CG warning: Matrix not PD");
			iterate.addInPlace(defect.mul(alpha));
			newResiduum = residuum.sub(z.mul(alpha));
			beta = newResiduum.inner(newResiduum) / residuum.inner(residuum);
			defect = newResiduum.add(defect.mul(beta));
			residuum = newResiduum;
			metric.publishIterate(residuum.euclidianNorm());
			if (showProgress)
			{
				System.out.println("CG residuum norm: " + residuum.euclidianNorm());
			}
		}
		i.running = false;
		ex.shutdown();
		lastSolution = iterate;
		return iterate;
	}
	
	public Vector solvePCG(final VectorMultiplyable operator,
	                       final VectorMultiplyable preconditioner,
	                       final Vector rhs,
	                       final double tol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			final DenseVector d = new DenseVector(rhs.getLength());
			for (int i = 0; i < d.getLength(); i++)
				d.set(Math.random(), i);
			if (operator.mvMul(d)
			            .absMaxElement() < 1e-14)
				System.err.println("operator has a kernel");
			if (preconditioner.mvMul(d)
			                  .absMaxElement() < 1e-14)
				System.err.println("preconditioner has a kernel");
		}
		metric.goal = tol;
		metric.restart();
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showInterrupt)
			ex.execute(i);
		final int n = rhs.getLength();
		Vector Ap;
		Vector z, newZ;
		Vector newResiduum;
		DenseVector iterate = new DenseVector(n);
		iterate.set(1, 0);
		if (lastSolution != null)
			iterate = new DenseVector(lastSolution);
		Vector residuum = rhs.sub(operator.mvMul(iterate));
		z = preconditioner.mvMul(residuum);
		Vector p = new DenseVector(z);
		double alpha;
		double beta;
		for (int iter = 0; iter < n && residuum.euclidianNorm() > tol && i.running; iter++)
		{
			Ap = operator.mvMul(p);
			alpha = residuum.inner(z) / p.inner(Ap);
			if (alpha < 0 && showProgress)
				System.out.println("CG warning: Matrix not PD");
			iterate.addInPlace(p.mul(alpha));
			newResiduum = residuum.sub(Ap.mul(alpha));
			newZ = preconditioner.mvMul(newResiduum);
			beta = newResiduum.inner(newZ) / residuum.inner(z);
			p = newZ.add(p.mul(beta));
			residuum = newResiduum;
			z = newZ;
			metric.publishIterate(residuum.euclidianNorm());
			if (showProgress)
				System.out.println("PCG residuum norm: " + residuum.euclidianNorm());
		}
		i.running = false;
		ex.shutdown();
		
		lastSolution = iterate;
		return iterate;
	}
	
	public <T extends VectorMultiplyable> Vector solvePGMRES(final VectorMultiplyable operator,
	                                                         final T preconditioner,
	                                                         final Vector rhs,
	                                                         final double tol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			final DenseVector d = new DenseVector(rhs.getLength());
			for (int i = 0; i < d.getLength(); i++)
				d.set(Math.random(), i);
			if (operator.mvMul(d)
			            .absMaxElement() < 1e-14)
				System.err.println("operator has a kernel");
			if (preconditioner.mvMul(d)
			                  .absMaxElement() < 1e-14)
				System.err.println("preconditioner has a kernel");
		}
		metric.goal = tol;
		metric.restart();
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showInterrupt)
			ex.execute(i);
		final Vector x;
		gm = new GMRES(showProgress, metric, 0);
		if (lastSolution == null)
			x = gm.solve(preconditioner, operator, rhs, tol, i);
		else
			x = gm.solve(preconditioner, operator, rhs, lastSolution, tol, i);
		i.running = false;
		ex.shutdown();
		lastSolution = x;
		return x;
	}
	
	public Vector solveGMRES(final VectorMultiplyable operator, final Vector rightHandSide, final double tol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			final DenseVector d = new DenseVector(rightHandSide.getLength());
			for (int i = 0; i < d.getLength(); i++)
				d.set(Math.random(), i);
			if (operator.mvMul(d)
			            .absMaxElement() < 1e-14)
				System.err.println("operator has a kernel");
		}
		metric.goal = tol;
		metric.restart();
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showInterrupt)
			ex.execute(i);
		final Vector x;
		gm = new GMRES(showProgress, metric, 0);
		if (lastSolution == null)
			x = gm.solve(operator, rightHandSide, tol, i);
		else
			x = gm.solve(operator, rightHandSide, lastSolution, tol, i);
		i.running = false;
		ex.shutdown();
		lastSolution = x;
		return x;
	}
	
	public Vector solveBiCGStab(final VectorMultiplyable operator, final Vector rhs, final double tol)
	{
		DenseVector iterate = new DenseVector(rhs.getLength());
		if (lastSolution != null)
			iterate = new DenseVector(lastSolution);
		else
			iterate.set(1, 0);
		return solveBiCGStab(operator, rhs, iterate, tol);
	}
	
	public Vector solvePBiCGStab(final VectorMultiplyable operator,
	                             final VectorMultiplyable preconditioner, final Vector rhs,
	                             final double tol)
	{
		DenseVector iterate = new DenseVector(rhs.getLength());
		if (lastSolution != null)
			iterate = new DenseVector(lastSolution);
		else
			iterate.set(1, 0);
		return solvePBiCGStab(operator, preconditioner, rhs, iterate, tol);
	}
	
	public Vector solveBiCGStab(final VectorMultiplyable operator,
	                            final Vector rhs,
	                            final Vector startIterate,
	                            final double tol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			final DenseVector d = new DenseVector(rhs.getLength());
			for (int i = 0; i < d.getLength(); i++)
				d.set(Math.random(), i);
			if (operator.mvMul(d)
			            .absMaxElement() < 1e-14)
				System.err.println("operator has a kernel");
		}
		metric.goal = tol;
		metric.restart();
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showInterrupt)
			ex.execute(i);
		final int n = rhs.getLength();
		Vector v = new DenseVector(n);
		Vector s;
		Vector t;
		DenseVector iterate = new DenseVector(startIterate);
		final Vector startResiduum = rhs.sub(operator.mvMul(iterate));
		Vector residuum = new DenseVector(startResiduum);
		Vector p = new DenseVector(n);
		double alpha = 1;
		double omega = 1;
		double rho = 1;
		double rhoLast = 1;
		double beta = 1;
		for (int iter = 0; iter < n && residuum.euclidianNorm() > tol && i.running; iter++)
		{
			rho = residuum.inner(startResiduum);
			if (Math.abs(rhoLast) < tol)
				return solveGMRES(operator, rhs, tol);
			beta = alpha / omega * rho / rhoLast;
			p = residuum.add(p.mul(beta))
			            .sub(v.mul(omega * beta));
			v = operator.mvMul(p);
			alpha = rho / v.inner(startResiduum);
			s = residuum.sub(v.mul(alpha));
			if (s.euclidianNorm() < tol)
			{
				if (showProgress)
					System.out.println("BiCGStab residuum norm: " + s.euclidianNorm());
				iterate = iterate.add(p.mul(alpha));
				break;
			}
			t = operator.mvMul(s);
			omega = t.inner(s) / t.inner(t);
			iterate = iterate.add(p.mul(alpha))
			                 .add(s.mul(omega));
			residuum = s.sub(t.mul(omega));
			rhoLast = rho;
			metric.publishIterate(residuum.euclidianNorm());
			if (showProgress)
				System.out.println("BiCGStab residuum norm: " + residuum.euclidianNorm());
		}
		i.running = false;
		ex.shutdown();
		lastSolution = iterate;
		return iterate;
	}
	
	public <T extends VectorMultiplyable> Vector solvePBiCGStab(final VectorMultiplyable operator,
	                                                            final T preconditioner,
	                                                            final Vector rhs,
	                                                            final Vector startIterate,
	                                                            final double tol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			final DenseVector d = new DenseVector(rhs.getLength());
			for (int i = 0; i < d.getLength(); i++)
				d.set(Math.random(), i);
			if (operator.mvMul(d)
			            .absMaxElement() < 1e-14)
				System.err.println("operator has a kernel");
			if (preconditioner.mvMul(d)
			                  .absMaxElement() < 1e-14)
				System.err.println("preconditioner has a kernel");
		}
		metric.goal = tol;
		metric.restart();
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showInterrupt)
			ex.execute(i);
		final int n = rhs.getLength();
		Vector s;
		Vector t;
		DenseVector iterate = new DenseVector(startIterate);
		final Vector startResiduum = rhs.sub(operator.mvMul(iterate));
		Vector residuum = new DenseVector(startResiduum);
		double alpha = 1;
		double omega = 1;
		double rho = 1;
		double rhoLast = 1;
		double beta = 1;
		Vector v = new DenseVector(n);
		Vector y = new DenseVector(n);
		Vector z = new DenseVector(n);
		Vector pt = new DenseVector(n);
		Vector p = new DenseVector(n);
		for (int iter = 0; iter < n * 10 && residuum.euclidianNorm() > tol && i.running; iter++)
		{
			rho = residuum.inner(startResiduum);
			if (Math.abs(rhoLast) < tol)
				return solvePGMRES(operator, preconditioner, rhs, tol);
			beta = alpha / omega * rho / rhoLast;
			p = residuum.add(p.mul(beta))
			            .sub(v.mul(omega * beta));
			y = preconditioner.mvMul(p);
			v = operator.mvMul(y);
			alpha = rho / v.inner(startResiduum);
			s = residuum.sub(v.mul(alpha));
			if (s.euclidianNorm() < tol)
			{
				if (showProgress)
					System.out.println("PBiCGStab residuum norm: " + s.euclidianNorm());
				iterate = iterate.add(y.mul(alpha));
				break;
			}
			z = preconditioner.mvMul(s);
			t = operator.mvMul(z);
			pt = preconditioner.mvMul(t);
			omega = pt.inner(z) / pt.inner(pt);
			iterate = iterate.add(y.mul(alpha))
			                 .add(z.mul(omega));
			residuum = s.sub(t.mul(omega));
			rhoLast = rho;
			metric.publishIterate(residuum.euclidianNorm());
			if (showProgress)
				System.out.println("PBiCGStab residuum norm: " + residuum.euclidianNorm());
		}
		i.running = false;
		ex.shutdown();
		lastSolution = iterate;
		return iterate;
	}
}
