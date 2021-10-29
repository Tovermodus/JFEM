package linalg;

import java.util.ArrayList;
import java.util.stream.IntStream;

class GMRES
{
	private final boolean showProgress;
	int restarts = 0;
	int MAX_RESTARTS = 100;
	int ITERATIONS_BEFORE_RESTART = 200;
	
	public GMRES(final boolean showProgress, final int restarts)
	{
		this.showProgress = showProgress;
		this.restarts = restarts;
	}
	
	public linalg.Vector solve(final VectorMultiplyable A,
	                           final Vector b,
	                           final double tol, final IterativeSolver.Interruptor interruptor)
	{
		
		final MutableVector x = new linalg.DenseVector(b.getLength());
		x.set(1, 0);
		return solve(A, b, x, tol, interruptor);
	}
	
	public linalg.Vector solve(final VectorMultiplyable A,
	                           final Vector b, final Vector x,
	                           final double tol, final IterativeSolver.Interruptor interruptor)
	{
		final ArrayList<Vector> v = new ArrayList<>();
		final int n = b.getLength();
		final MutableVector c = new linalg.DenseVector(n);
		final MutableVector s = new linalg.DenseVector(n);
		final MutableVector gamma = new linalg.DenseVector(n + 1);
		final Vector r = b.sub(A.mvMul(x));
		final DenseMatrix h = new DenseMatrix(110, 110);
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
					 .mapToDouble(vec -> vec.inner(q))
					 .toArray());
			h.addSmallMatrixInPlaceAt(newHValues.asMatrix(), 0, j);
			final Vector newHV =
				(IntStream.range(0, j + 1))
					.mapToObj(i -> v.get(i)
					                .mul(newHValues.at(i)))
					.reduce(new DenseVector(n), Vector::add);
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
			if (showProgress)
				System.out.println(Math.abs(gamma.at(j + 1)));
			if (Math.abs(gamma.at(j + 1)) < tol || j > ITERATIONS_BEFORE_RESTART)
				break;
			v.add(w.mul(1. / h.at(j + 1, j)));
		}
		if (j == n || !interruptor.isRunning())
			j--;
		final DenseVector alpha = new DenseVector(n);
		for (int i = j; i >= 0; i--)
		{
			final int finalI = i;
			final double hAlpha = IntStream
				.range(i + 1, j + 1)
				.mapToDouble(k -> h.at(finalI, k) * alpha.at(k))
				.sum();
			alpha.set(1. / h.at(i, i) * (gamma.at(i) - hAlpha), i);
		}
		final Vector alphaV =
			IntStream.range(0, j + 1)
			         .mapToObj(i -> v.get(i)
			                         .mul(alpha.at(i)))
			         .reduce(
				         new DenseVector(n),
				         Vector::add);
		if (j > ITERATIONS_BEFORE_RESTART && restarts < MAX_RESTARTS)
		{
			restarts++;
			return solve(A, b, x.add(alphaV), tol, interruptor);
		}
		return x.add(alphaV);
	}
	
	public <T extends VectorMultiplyable> linalg.Vector solve(final T preconditioner,
	                                                          final VectorMultiplyable A,
	                                                          final Vector b,
	                                                          final double tol,
	                                                          final IterativeSolver.Interruptor interruptor)
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
	                                                          final IterativeSolver.Interruptor interruptor)
	{
		final ArrayList<Vector> v = new ArrayList<>();
		final int n = b.getLength();
		final MutableVector c = new linalg.DenseVector(n);
		final MutableVector s = new linalg.DenseVector(n);
		final MutableVector gamma = new linalg.DenseVector(n + 1);
		final Vector r = preconditioner.mvMul(b.sub(A.mvMul(x)));
		final SparseMatrix h = new SparseMatrix(n + 1, n + 1);
		if (b.sub(A.mvMul(x))
		     .euclidianNorm() <= tol)
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
					 .mapToDouble(vec -> vec.inner(q))
					 .toArray());
			h.addSmallMatrixAt(newHValues.asMatrix(), 0, j);
			final Vector newHV =
				(IntStream.range(0, j + 1))
					.mapToObj(i -> v.get(i)
					                .mul(newHValues.at(i)))
					.reduce(new DenseVector(n), Vector::add);
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
			if (showProgress)
				System.out.println(Math.abs(gamma.at(j + 1)));
			if (Math.abs(gamma.at(j + 1)) < tol || j > ITERATIONS_BEFORE_RESTART)
				break;
			v.add(w.mul(1. / h.at(j + 1, j)));
		}
		if (j == n || !interruptor.isRunning())
			j--;
		final DenseVector alpha = new DenseVector(n);
		for (int i = j; i >= 0; i--)
		{
			final int finalI = i;
			final double hAlpha = IntStream
				.range(i + 1, j + 1)
				.mapToDouble(k -> h.at(finalI, k) * alpha.at(k))
				.sum();
			alpha.set(1. / h.at(i, i) * (gamma.at(i) - hAlpha), i);
		}
		final Vector alphaV =
			IntStream.range(0, j + 1)
			         .mapToObj(i -> v.get(i)
			                         .mul(alpha.at(i)))
			         .reduce(
				         new DenseVector(n),
				         Vector::add);
		if (j > ITERATIONS_BEFORE_RESTART && restarts < MAX_RESTARTS)
		{
			restarts++;
			return solve(preconditioner, A, b, x.add(alphaV), tol, interruptor);
		}
		return x.add(alphaV);
	}
}
