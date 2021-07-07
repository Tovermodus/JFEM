package linalg;

import com.google.common.base.Stopwatch;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class IterativeSolver<Op extends VectorMultiplyable>
{
	private volatile boolean interrupted = false;
	private GMRES gm;
	ExecutorService ex;
	public boolean showProgress = true;
	public Vector solveCG(Op operator, Vector rhs, double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		Interruptor i = new Interruptor();
		if(showProgress)
			ex.execute(i);
		Stopwatch s = Stopwatch.createStarted();
		int n = rhs.getLength();
		Vector z;
		Vector newResiduum;
		DenseVector iterate = new DenseVector(n);
		iterate.set(Math.random(),0);
		Vector m = operator.mvMul(iterate);
		Vector residuum = rhs.sub(m);
		Vector defect = new DenseVector(residuum);
		double alpha;
		double beta;
		for(int iter = 0; iter < n && residuum.euclidianNorm() > tol && !interrupted; iter++)
		{
			z = operator.mvMul(defect);
			alpha = residuum.inner(residuum)/defect.inner(z);
			iterate.addInPlace(defect.mul(alpha));
			newResiduum = residuum.sub(z.mul(alpha));
			beta = newResiduum.inner(newResiduum)/residuum.inner(residuum);
			defect = newResiduum.add(defect.mul(beta));
			residuum = newResiduum;
			if(showProgress)
			System.out.println(residuum.euclidianNorm());
		}
		interrupted = true;
		ex.shutdown();
		
		return iterate;
	}
	public <T extends VectorMultiplyable> Vector solvePGMRES(Op operator, T preconditioner, Vector rhs,
	                                                              double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		Interruptor i = new Interruptor();
		if(showProgress)
		ex.execute(i);
		int n = rhs.getLength();
		DenseVector v = new DenseVector(n);
		Vector x = null;
		gm = new GMRES(v.toMTJvector());
		try
		{
			x =  gm.solve(preconditioner, operator, rhs, tol);
		} catch (IterativeSolverNotConvergedException e)
		{
			e.printStackTrace();
		}
		interrupted = true;
		ex.shutdown();
		return x;
	}
	public Vector solveGMRES(Op operator, Vector rightHandSide, double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		Interruptor i = new Interruptor();
		if(showProgress)
		ex.execute(i);
		int n = rightHandSide.getLength();
		DenseVector v = new DenseVector(n);
		Vector x = null;
		gm = new GMRES(v.toMTJvector());
		try
		{
			x =  gm.solve(operator, rightHandSide, tol);
		} catch (IterativeSolverNotConvergedException e)
		{
			e.printStackTrace();
		}
		interrupted = true;
		ex.shutdown();
		return x;
	}
	public Vector solveBiCGStab(Op operator, Vector rhs, double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		Interruptor i = new Interruptor();
		if(showProgress)
		ex.execute(i);
		int n = rhs.getLength();
		Vector v;
		Vector s;
		Vector t;
		Vector iterate = new DenseVector(n);
		Vector startResiduum = rhs.sub(operator.mvMul(iterate));
		Vector residuum = new DenseVector(startResiduum);
		Vector p = new DenseVector(residuum);
		double alpha;
		double omega;
		double rho = residuum.inner(residuum);
		double rhoLast;
		double beta;
		for(int iter = 0; iter < n && residuum.euclidianNorm() > tol&& !interrupted; iter++)
		{
			v = operator.mvMul(p);
			alpha = rho/v.inner(startResiduum);
			s = residuum.sub(v.mul(alpha));
			t = operator.mvMul(s);
			omega = t.inner(s)/t.inner(t);
			iterate = iterate.add(p.mul(alpha)).add(s.mul(omega));
			residuum = s.sub(t.mul(omega));
			rhoLast = rho;
			rho = residuum.inner(startResiduum);
			beta = alpha/omega*rho/rhoLast;
			p = residuum.add(p.mul(beta)).sub(v.mul(omega*beta));
			if(showProgress)
			System.out.println(residuum.euclidianNorm());
		}
		interrupted = true;
		ex.shutdown();
		return iterate;


	}
	public <T extends VectorMultiplyable> Vector solvePBiCGStab(Op operator, T preconditioner, Vector rhs,
	                                                                 double tol)
	{
		ex = Executors.newSingleThreadExecutor();
		Interruptor i = new Interruptor();
		if(showProgress)
		ex.execute(i);
		int n = rhs.getLength();
		Vector v;
		Vector vP;
		Vector s;
		Vector sP;
		Vector t;
		Vector tP;
		Vector iterate = new DenseVector(n);
		Vector startResiduum = rhs.sub(operator.mvMul(iterate));
		Vector startResiduumP = preconditioner.mvMul(startResiduum);
		Vector residuum = new DenseVector(startResiduum);
		Vector residuumP = preconditioner.mvMul(residuum);
		Vector pP = preconditioner.mvMul(residuumP);
		double alpha;
		double omega;
		double rho = residuumP.inner(residuumP);
		double rhoLast;
		double beta;
		for(int iter = 0; iter < n && residuum.euclidianNorm() > tol&& !interrupted; iter++)
		{
			v = operator.mvMul(pP);
			vP = preconditioner.mvMul(v);
			alpha = rho/vP.inner(startResiduumP);
			s = residuum.sub(v.mul(alpha));
			sP = preconditioner.mvMul(s);
			t = operator.mvMul(sP);
			tP = preconditioner.mvMul(t);
			omega = tP.inner(sP)/tP.inner(tP);
			iterate = iterate.add(pP.mul(alpha)).add(sP.mul(omega));
			residuum = rhs.sub(operator.mvMul(iterate));
			residuumP = sP.sub(tP.mul(omega));
			rhoLast = rho;
			rho = residuumP.inner(startResiduumP);
			beta = alpha/omega*rho/rhoLast;
			pP = residuumP.add(pP.mul(beta)).sub(vP.mul(omega*beta));
		}
		interrupted = true;
		ex.shutdown();
		return iterate;


	}
	class Interruptor implements Runnable
	{
		private JFrame f;
		@Override
		public void run()
		{
			f = new JFrame("Interrupt Solver");
			JButton b = new JButton("Interrupt!");
			f.setLayout(new GridLayout());
			f.add(b);
			f.setBounds(0,0,300,300);
			f.setVisible(true);
			b.addActionListener(new ActionListener()
			{
				@Override
				public void actionPerformed(ActionEvent e)
				{
					System.out.println("Interrupt!!!");
					interrupted = true;
					if(gm != null)
						gm.interrupted = true;
					f.setVisible(false);
					f.dispose();
				}
			});
			while(true)
			{
				if(interrupted)
					f.dispose();
			}
		}
		
	}
}
