package linalg;

import basic.Interruptor;
import basic.PerformanceArguments;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.function.UnaryOperator;
import java.util.stream.IntStream;

public class GMRES
{
	public boolean printProgress;
	final IterativeSolverConvergenceMetric metric;
	public int MAX_RESTARTS = 200;
	int restarts;
	int iterations;
	public int ITERATIONS_BEFORE_RESTART = 150;
	
	public GMRES(final boolean printProgress,
	             final IterativeSolverConvergenceMetric metric)
	{
		this.printProgress = printProgress;
		this.restarts = 0;
		this.metric = metric;
		iterations = 0;
	}
	
	public GMRES(final GMRES old)
	{
		this.printProgress = old.printProgress;
		this.restarts = old.restarts + 1;
		this.metric = old.metric;
		this.MAX_RESTARTS = old.MAX_RESTARTS;
		this.ITERATIONS_BEFORE_RESTART = old.ITERATIONS_BEFORE_RESTART;
		iterations = 0;
	}
	
	public linalg.Vector solve(final VectorMultiplyable A,
	                           final Vector b,
	                           final double tol, final Interruptor interruptor)
	{
		
		final MutableVector x = new DenseVector(b.getLength());
		x.set(1, 0);
		return solve(A, b, x, tol, interruptor);
	}
	
	public Vector solve(final VectorMultiplyable preconditioner,
	                    final VectorMultiplyable A,
	                    final Vector b,
	                    final double tol,
	                    final Interruptor interruptor)
	{
		
		final MutableVector x = new DenseVector(b.getLength());
		x.set(1, 0);
		return solve(preconditioner, A, b, x, tol, interruptor);
	}
	
	public linalg.Vector solve(final VectorMultiplyable preconditioner,
	                           final VectorMultiplyable A,
	                           final Vector b,
	                           final Vector x,
	                           final double tol,
	                           final Interruptor interruptor)
	{
		metric.goal = tol;
		final ArrayList<Vector> v = new ArrayList<>();
		final int n = b.getLength();
		final MutableVector gamma = new linalg.DenseVector(ITERATIONS_BEFORE_RESTART + 10);
		final Vector r = preconditioner.mvMul(b.sub(A.mvMul(x)));
		final SparseMatrix h = new SparseMatrix(ITERATIONS_BEFORE_RESTART + 10, ITERATIONS_BEFORE_RESTART + 10);
		if (r.euclidianNorm() <= tol)
			return x;
		gamma.set(r.euclidianNorm(), 0);
		v.add(r.mul(1. / gamma.at(0)));
		int j = innerGMRES(vec -> preconditioner.mvMul(A.mvMul(vec)), A, b, tol, v, gamma, h, x,
		                   n, interruptor);
		if (j == n || !interruptor.isRunning())
			j--;
		final DenseVector alpha = calculateAlpha(gamma, h, j);
		final Vector alphaV =
			linearCombination(v, alpha);
		if (j > ITERATIONS_BEFORE_RESTART && restarts < MAX_RESTARTS)
		{
			restarts++;
			if (printProgress)
				System.out.println("GMRES restart");
			return solve(preconditioner, A, b, x.add(alphaV), tol, interruptor);
		}
		return x.add(alphaV);
	}
	
	public linalg.Vector solve(final VectorMultiplyable A,
	                           final Vector b, final Vector x,
	                           final double tol, final Interruptor interruptor)
	{
		metric.goal = tol;
		final ArrayList<Vector> v = new ArrayList<>();
		final int n = b.getLength();
		final MutableVector gamma = new linalg.DenseVector(ITERATIONS_BEFORE_RESTART + 10);
		final Vector r = b.sub(A.mvMul(x));
		final SparseMatrix h = new SparseMatrix(ITERATIONS_BEFORE_RESTART + 10, ITERATIONS_BEFORE_RESTART + 10);
		if (r.euclidianNorm() <= tol)
			return x;
		gamma.set(r.euclidianNorm(), 0);
		v.add(r.mul(1. / gamma.at(0)));
		int j = innerGMRES(A::mvMul, A, b, tol, v, gamma, h, x, n, interruptor);
		if (j == n || !interruptor.isRunning())
			j--;
		final DenseVector alpha = calculateAlpha(gamma, h, j);
		final Vector alphaV = linearCombination(v, alpha);
		if (j > ITERATIONS_BEFORE_RESTART && restarts < MAX_RESTARTS)
		{
			restarts++;
			if (printProgress)
				System.out.println("GMRes Restart");
			return solve(A, b, x.add(alphaV), tol, interruptor);
		}
		return x.add(alphaV);
	}
	
	private int innerGMRES(final UnaryOperator<Vector> vectorMultiplocation,
	                       final VectorMultiplyable A,
	                       final Vector b,
	                       final double tol,
	                       final ArrayList<Vector> v,
	                       final MutableVector gamma,
	                       final SparseMatrix h,
	                       final Vector x,
	                       final int n,
	                       final Interruptor interruptor)
	{
		final MutableVector c = new linalg.DenseVector(n);
		final MutableVector s = new linalg.DenseVector(n);
		int j;
		for (j = 0; j < n && interruptor.isRunning(); j++)
		{
			final Vector q = vectorMultiplocation.apply(v.get(j));
			final DenseVector orthogonalityFactors = calculateOrthogonalityFactors(v, q);
			final Vector orthogonalPartInSubspace = linearCombination(v, orthogonalityFactors);
			final Vector w = q.sub(orthogonalPartInSubspace);
			h.addColumn(orthogonalityFactors, j);
			h.add(w.euclidianNorm(), j + 1, j);
			calculateGivens(c, s, h, j);
			final double beta = Math.sqrt(Math.pow(h.at(j, j), 2) + Math.pow(h.at(j + 1, j), 2));
			s.set(h.at(j + 1, j) / beta, j);
			c.set(h.at(j, j) / beta, j);
			h.set(beta, j, j);
			gamma.set(-s.at(j) * gamma.at(j), j + 1);
			gamma.set(c.at(j) * gamma.at(j), j);
			if (printProgress)
				System.out.println("GMRes Gamma: " + Math.abs(gamma.at(j + 1)));
			saveResidual(A, b, v, gamma, h, x, j);
			if (Math.abs(gamma.at(j + 1)) < tol || j > ITERATIONS_BEFORE_RESTART)
				return j;
			v.add(w.mul(1. / h.at(j + 1, j)));
			iterations++;
		}
		return j;
	}
	
	private void saveResidual(final VectorMultiplyable A,
	                          final Vector b,
	                          final ArrayList<Vector> v,
	                          final MutableVector gamma,
	                          final SparseMatrix h,
	                          final Vector x,
	                          final int j)
	{
		if (PerformanceArguments.getInstance().GMResData == PerformanceArguments.GMRESResidual)
			metric.publishIterateAsync(() ->
			                           {
				                           final DenseVector alpha = calculateAlpha(gamma, h, j);
				                           final Vector step = linearCombination(v, alpha);
				                           final Vector newIt = x.add(step);
				                           return A.mvMul(newIt)
				                                   .sub(b)
				                                   .euclidianNorm();
			                           });
		else
			metric.publishIterate(gamma.at(j + 1));
	}
	
	@NotNull
	private static DenseVector calculateOrthogonalityFactors(final ArrayList<Vector> v, final Vector q)
	{
		return DenseVector.vectorFromValues(
			v.stream()
			 .parallel()
			 .mapToDouble(vec -> vec.inner(q))
			 .toArray());
	}
	
	private static void calculateGivens(final MutableVector c,
	                                    final MutableVector s,
	                                    final SparseMatrix h,
	                                    final int j)
	{
		for (int i = 0; i < j; i++)
		{
			final double c_i = c.at(i);
			final double s_i = s.at(i);
			final double h_ij = h.at(i, j);
			final double h_ipj = h.at(i + 1, j);
			h.set(h_ij * c_i + h_ipj * s_i, i, j);
			h.set(-h_ij * s_i + h_ipj * c_i, i + 1, j);
		}
	}
	
	private static Vector linearCombination(final ArrayList<Vector> vectors,
	                                        final DenseVector coefficients)
	{
		final int n = vectors.get(0)
		                     .getLength();
		return IntStream.range(0, coefficients.getLength())
		                .mapToObj(i -> vectors.get(i)
		                                      .mul(coefficients.at(i)))
		                .reduce(new DenseVector(n), Vector::add);
	}
	
	private static DenseVector calculateAlpha(final MutableVector gamma, final Matrix h, final int j)
	{
		final DenseVector alpha = new DenseVector(j + 1);
		for (int i = j; i >= 0; i--)
		{
			final int finalI = i;
			final double hAlpha = IntStream.range(i + 1, j + 1)
			                               .mapToDouble(k -> h.at(finalI, k) * alpha.at(k))
			                               .sum();
			alpha.set(1. / h.at(i, i) * (gamma.at(i) - hAlpha), i);
		}
		return alpha;
	}
}
