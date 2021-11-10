package linalg;

import basic.Interruptor;

import java.util.ArrayList;
import java.util.stream.IntStream;

public class GMRES
{
	private final boolean printProgress;
	final IterativeSolverConvergenceMetric metric;
	public int MAX_RESTARTS = 200;
	int restarts = 0;
	public int ITERATIONS_BEFORE_RESTART = 150;
	
	public GMRES(final boolean printProgress,
	             final IterativeSolverConvergenceMetric metric, final int restarts)
	{
		this.printProgress = printProgress;
		this.restarts = restarts;
		this.metric = metric;
	}
	
	public linalg.Vector solve(final VectorMultiplyable A,
	                           final Vector b,
	                           final double tol, final Interruptor interruptor)
	{
		
		final MutableVector x = new linalg.DenseVector(b.getLength());
		x.set(1, 0);
		return solve(A, b, x, tol, interruptor);
	}
	
	public linalg.Vector solve(final VectorMultiplyable A,
	                           final Vector b, final Vector x,
	                           final double tol, final Interruptor interruptor)
	{
		metric.goal = tol;
		final ArrayList<Vector> v = new ArrayList<>();
		final int n = b.getLength();
		final MutableVector c = new linalg.DenseVector(n);
		final MutableVector s = new linalg.DenseVector(n);
		final MutableVector gamma = new linalg.DenseVector(n + 1);
		final Vector r = b.sub(A.mvMul(x));
		final DenseMatrix h = new DenseMatrix(ITERATIONS_BEFORE_RESTART + 10, ITERATIONS_BEFORE_RESTART + 10);
		if (r.euclidianNorm() <= tol)
			return x;
		gamma.set(r.euclidianNorm(), 0);
		v.add(r.mul(1. / gamma.at(0)));
		int j;
		for (j = 0; j < n && r.euclidianNorm() > tol && interruptor.isRunning(); j++)
		{
			final Vector q = A.mvMul(v.get(j));
			final DenseVector newHValues =
				DenseVector.vectorFromValues(
					v.stream()
					 .parallel()
					 .mapToDouble(vec -> vec.inner(q))
					 .toArray());
			h.addSmallMatrixInPlaceAt(newHValues.asMatrix(), 0, j);
			final Vector newHV =
				calculateStep(v, j, newHValues);
			final Vector w = q.sub(newHV);
			h.add(w.euclidianNorm(), j + 1, j);
			for (int i = 0; i < j; i++)
			{
				final double c_i = c.at(i);
				final double s_i = s.at(i);
				final double h_ij = h.at(i, j);
				final double h_ipj = h.at(i + 1, j);
				h.set(h_ij * c_i + h_ipj * s_i, i, j);
				h.set(-h_ij * s_i + h_ipj * c_i, i + 1, j);
			}
			final double beta = Math.sqrt(Math.pow(h.at(j, j), 2) + Math.pow(h.at(j + 1, j), 2));
			s.set(h.at(j + 1, j) / beta, j);
			c.set(h.at(j, j) / beta, j);
			h.set(beta, j, j);
			gamma.set(-s.at(j) * gamma.at(j), j + 1);
			gamma.set(c.at(j) * gamma.at(j), j);
			if (printProgress)
				System.out.println("GMRes Gamma: " + Math.abs(gamma.at(j + 1)));
			final int finalJ = j;
			metric.publishIterateAsync(() ->
			                           {
				                           final DenseVector alpha = calculateAlpha(gamma, h, finalJ);
				                           final Vector step = calculateStep(v, finalJ, alpha);
				                           final Vector newIt = x.add(step);
				                           return A.mvMul(newIt)
				                                   .sub(b)
				                                   .euclidianNorm();
			                           });
			if (Math.abs(gamma.at(j + 1)) < tol || j > ITERATIONS_BEFORE_RESTART)
				break;
			v.add(w.mul(1. / h.at(j + 1, j)));
		}
		if (j == n || !interruptor.isRunning())
			j--;
		final DenseVector alpha = calculateAlpha(gamma, h, j);
		final Vector alphaV =
			calculateStep(v, j, alpha);
		if (j > ITERATIONS_BEFORE_RESTART && restarts < MAX_RESTARTS)
		{
			restarts++;
			if (printProgress)
				System.out.println("GMRes Restart");
			return solve(A, b, x.add(alphaV), tol, interruptor);
		}
		return x.add(alphaV);
	}
	
	public <T extends VectorMultiplyable> linalg.Vector solve(final T preconditioner,
	                                                          final VectorMultiplyable A,
	                                                          final Vector b,
	                                                          final double tol,
	                                                          final Interruptor interruptor)
	{
		
		final MutableVector x = new linalg.DenseVector(b.getLength());
		x.set(1, 0);
		return solve(preconditioner, A, b, x, tol, interruptor);
	}
	
	public <T extends VectorMultiplyable> linalg.Vector solve(final T preconditioner,
	                                                          final VectorMultiplyable A,
	                                                          final Vector b,
	                                                          final Vector x,
	                                                          final double tol,
	                                                          final Interruptor interruptor)
	{
		metric.goal = tol;
		final ArrayList<Vector> v = new ArrayList<>();
		final int n = b.getLength();
		final MutableVector c = new linalg.DenseVector(n);
		final MutableVector s = new linalg.DenseVector(n);
		final MutableVector gamma = new linalg.DenseVector(n + 1);
		final Vector r = preconditioner.mvMul(b.sub(A.mvMul(x)));
		final SparseMatrix h = new SparseMatrix(ITERATIONS_BEFORE_RESTART + 10, ITERATIONS_BEFORE_RESTART + 10);
		if (r.euclidianNorm() <= tol)
			return x;
		gamma.set(r.euclidianNorm(), 0);
		v.add(r.mul(1. / gamma.at(0)));
		int j;
		for (j = 0; j < n && r.euclidianNorm() > tol && interruptor.isRunning(); j++)
		{
			final Vector q = preconditioner.mvMul(A.mvMul(v.get(j)));
			final DenseVector newHValues =
				DenseVector.vectorFromValues(
					v.stream()
					 .parallel()
					 .mapToDouble(vec -> vec.inner(q))
					 .toArray());
			h.addSmallMatrixAt(newHValues.asMatrix(), 0, j);
			final Vector newHV =
				calculateStep(v, j, newHValues);
			final Vector w = q.sub(newHV);
			h.add(w.euclidianNorm(), j + 1, j);
			for (int i = 0; i < j; i++)
			{
				final double c_i = c.at(i);
				final double s_i = s.at(i);
				final double h_ij = h.at(i, j);
				final double h_ipj = h.at(i + 1, j);
				h.set(h_ij * c_i + h_ipj * s_i, i, j);
				h.set(-h_ij * s_i + h_ipj * c_i, i + 1, j);
			}
			final double beta = Math.sqrt(Math.pow(h.at(j, j), 2) + Math.pow(h.at(j + 1, j), 2));
			s.set(h.at(j + 1, j) / beta, j);
			c.set(h.at(j, j) / beta, j);
			h.set(beta, j, j);
			gamma.set(-s.at(j) * gamma.at(j), j + 1);
			gamma.set(c.at(j) * gamma.at(j), j);
			if (printProgress)
				System.out.println("GMRes Gamma: " + Math.abs(gamma.at(j + 1)));
			//metric.publishIterate(Math.abs(gamma.at(j + 1)));
			final int finalJ = j;
			metric.publishIterateAsync(() ->
			                           {
				                           final DenseVector alpha = calculateAlpha(gamma, h, finalJ);
				                           final Vector step = calculateStep(v, finalJ, alpha);
				                           final Vector newIt = x.add(step);
				                           return A.mvMul(newIt)
				                                   .sub(b)
				                                   .euclidianNorm();
			                           });
			if (Math.abs(gamma.at(j + 1)) < tol || j > ITERATIONS_BEFORE_RESTART)
				break;
			v.add(w.mul(1. / h.at(j + 1, j)));
		}
		if (j == n || !interruptor.isRunning())
			j--;
		final DenseVector alpha = calculateAlpha(gamma, h, j);
		final Vector alphaV =
			calculateStep(v, j, alpha);
		if (j > ITERATIONS_BEFORE_RESTART && restarts < MAX_RESTARTS)
		{
			restarts++;
			if (printProgress)
				System.out.println("GMRES restart");
			return solve(preconditioner, A, b, x.add(alphaV), tol, interruptor);
		}
		return x.add(alphaV);
	}
	
	private static Vector calculateStep(final ArrayList<Vector> v,
	                                    final int j,
	                                    final DenseVector alpha)
	{
		final int n = v.get(0)
		               .getLength();
		return IntStream.range(0, j + 1)
		                .mapToObj(i -> v.get(i)
		                                .mul(alpha.at(i)))
		                .reduce(new DenseVector(n), Vector::add);
	}
	
	private static DenseVector calculateAlpha(final MutableVector gamma, final Matrix h, final int j)
	{
		final DenseVector alpha = new DenseVector(j + 1);
		for (int i = j; i >= 0; i--)
		{
			final int finalI = i;
			final double hAlpha = IntStream
				.range(i + 1, j + 1)
				.mapToDouble(k -> h.at(finalI, k) * alpha.at(k))
				.sum();
			alpha.set(1. / h.at(i, i) * (gamma.at(i) - hAlpha), i);
		}
		return alpha;
	}
}
