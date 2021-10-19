package linalg;

import javax.swing.*;
import java.awt.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class IterativeSolver
{
	private GMRES gm;
	ExecutorService ex;
	public boolean showProgress = true;
	
	public Vector solveCG(final VectorMultiplyable operator, final Vector rhs, final double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showProgress)
			ex.execute(i);
		final int n = rhs.getLength();
		Vector z;
		Vector newResiduum;
		final DenseVector iterate = new DenseVector(n);
		iterate.set(1, 0);
		final Vector m = operator.mvMul(iterate);
		Vector residuum = rhs.sub(m);
		Vector defect = new DenseVector(residuum);
		double alpha;
		double beta;
		for (int iter = 0; iter < n && residuum.euclidianNorm() > tol && i.running; iter++)
		{
			z = operator.mvMul(defect);
			alpha = residuum.inner(residuum) / defect.inner(z);
			iterate.addInPlace(defect.mul(alpha));
			newResiduum = residuum.sub(z.mul(alpha));
			beta = newResiduum.inner(newResiduum) / residuum.inner(residuum);
			defect = newResiduum.add(defect.mul(beta));
			residuum = newResiduum;
			if (showProgress)
				System.out.println(residuum.euclidianNorm());
		}
		i.running = false;
		ex.shutdown();
		
		return iterate;
	}
	
	public Vector solvePCG(final VectorMultiplyable operator, final VectorMultiplyable preconditioner, final Vector rhs, final double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showProgress)
			ex.execute(i);
		final int n = rhs.getLength();
		Vector Ap;
		Vector z, newZ;
		Vector newResiduum;
		final DenseVector iterate = new DenseVector(n);
		iterate.set(1, 0);
		Vector residuum = rhs.sub(operator.mvMul(iterate));
		z = preconditioner.mvMul(residuum);
		Vector p = new DenseVector(z);
		double alpha;
		double beta;
		for (int iter = 0; iter < n && residuum.euclidianNorm() > tol && i.running; iter++)
		{
			Ap = operator.mvMul(p);
			alpha = residuum.inner(z) / p.inner(Ap);
			iterate.addInPlace(p.mul(alpha));
			newResiduum = residuum.sub(Ap.mul(alpha));
			newZ = preconditioner.mvMul(newResiduum);
			beta = newResiduum.inner(newZ) / residuum.inner(z);
			p = newZ.add(p.mul(beta));
			residuum = newResiduum;
			z = newZ;
			if (showProgress)
				System.out.println(residuum.euclidianNorm());
		}
		i.running = false;
		ex.shutdown();
		
		return iterate;
	}
	
	public <T extends VectorMultiplyable> Vector solvePGMRES(final VectorMultiplyable operator, final T preconditioner, final Vector rhs,
	                                                         final double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showProgress)
			ex.execute(i);
		final Vector x;
		gm = new GMRES(showProgress, 0);
		x = gm.solve(preconditioner, operator, rhs, tol, i);
		i.running = false;
		ex.shutdown();
		return x;
	}
	
	public Vector solveGMRES(final VectorMultiplyable operator, final Vector rightHandSide, final double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showProgress)
			ex.execute(i);
		final Vector x;
		gm = new GMRES(showProgress, 0);
		x = gm.solve(operator, rightHandSide, tol, i);
		i.running = false;
		ex.shutdown();
		return x;
	}
	
	public Vector solveBiCGStab(final VectorMultiplyable operator, final Vector rhs, final double tol)
	{
		return solveBiCGStab(operator, rhs, new DenseVector(rhs.getLength()), tol);
	}
	
	public Vector solveBiCGStab(final VectorMultiplyable operator, final Vector rhs, final Vector startIterate, final double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showProgress)
			ex.execute(i);
		final int n = rhs.getLength();
		Vector v;
		Vector s;
		Vector t;
		Vector iterate = new DenseVector(startIterate);
		final Vector startResiduum = rhs.sub(operator.mvMul(iterate));
		Vector residuum = new DenseVector(startResiduum);
		Vector p = new DenseVector(residuum);
		double alpha;
		double omega;
		double rho = residuum.inner(residuum);
		double rhoLast;
		double beta;
		for (int iter = 0; iter < n && residuum.euclidianNorm() > tol && i.running; iter++)
		{
			v = operator.mvMul(p);
			alpha = rho / v.inner(startResiduum);
			s = residuum.sub(v.mul(alpha));
			t = operator.mvMul(s);
			omega = t.inner(s) / t.inner(t);
			iterate = iterate.add(p.mul(alpha)).add(s.mul(omega));
			residuum = s.sub(t.mul(omega));
			rhoLast = rho;
			rho = residuum.inner(startResiduum);
			beta = alpha / omega * rho / rhoLast;
			p = residuum.add(p.mul(beta)).sub(v.mul(omega * beta));
			if (showProgress)
				System.out.println(residuum.euclidianNorm());
		}
		i.running = false;
		ex.shutdown();
		return iterate;
	}
	
	public <T extends VectorMultiplyable> Vector solvePBiCGStab(final VectorMultiplyable operator, final T preconditioner, final Vector rhs,
	                                                            final double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		final Interruptor i = new Interruptor();
		if (showProgress)
			ex.execute(i);
		final int n = rhs.getLength();
		Vector v;
		Vector vP;
		Vector s;
		Vector sP;
		Vector t;
		Vector tP;
		Vector iterate = new DenseVector(n);
		final Vector startResiduum = rhs.sub(operator.mvMul(iterate));
		final Vector startResiduumP = preconditioner.mvMul(startResiduum);
		Vector residuum = new DenseVector(startResiduum);
		Vector residuumP = preconditioner.mvMul(residuum);
		Vector pP = preconditioner.mvMul(residuumP);
		double alpha;
		double omega;
		double rho = residuumP.inner(residuumP);
		double rhoLast;
		double beta;
		for (int iter = 0; iter < n && residuum.euclidianNorm() > tol && i.running; iter++)
		{
			v = operator.mvMul(pP);
			vP = preconditioner.mvMul(v);
			alpha = rho / vP.inner(startResiduumP);
			s = residuum.sub(v.mul(alpha));
			sP = preconditioner.mvMul(s);
			t = operator.mvMul(sP);
			tP = preconditioner.mvMul(t);
			omega = tP.inner(sP) / tP.inner(tP);
			iterate = iterate.add(pP.mul(alpha)).add(sP.mul(omega));
			residuum = rhs.sub(operator.mvMul(iterate));
			residuumP = sP.sub(tP.mul(omega));
			rhoLast = rho;
			rho = residuumP.inner(startResiduumP);
			beta = alpha / omega * rho / rhoLast;
			pP = residuumP.add(pP.mul(beta)).sub(vP.mul(omega * beta));
		}
		i.running = false;
		ex.shutdown();
		return iterate;
	}
	
	public static class Interruptor implements Runnable
	{
		private JFrame f;
		public volatile boolean running = true;
		
		public boolean isRunning()
		{
			return running;
		}
		
		@Override
		public void run()
		{
			f = new JFrame("Interrupt Solver");
			final JButton b = new JButton("Interrupt!");
			f.setLayout(new GridLayout());
			f.add(b);
			f.setBounds(0, 0, 300, 300);
			f.setVisible(true);
			b.addActionListener(e ->
			                    {
				                    System.out.println("Interrupt!!!");
				                    running = false;
				                    f.setVisible(false);
				                    f.dispose();
			                    });
			while (running)
			{
				try
				{
					Thread.sleep(30);
				} catch (final InterruptedException e)
				{
					e.printStackTrace();
				}
			}
			f.dispose();
		}
	}
}
