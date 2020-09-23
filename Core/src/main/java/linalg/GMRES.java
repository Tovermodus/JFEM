package linalg;/*
 * Copyright (C) 2003-2006 Bjørn-Ove Heimsund
 *
 * This file is part of MTJ.
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
/*
 * Derived from public domain software at http://www.netlib.org/templates
 */


import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.GivensRotation;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.UpperTriangDenseMatrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.AbstractIterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;

/**
 * GMRES solver. GMRES solves the unsymmetric linear system <code>Ax = b</code>
 * using the Generalized Minimum Residual method. The GMRES iteration is
 * restarted after a given number of iterations. By default it is restarted
 * after 30 iterations.
 *
 * @author Templates
 */
class GMRES extends AbstractIterativeSolver
{

	/**
	 * After this many iterations, the GMRES will be restarted.
	 */
	private int restart;

	/**
	 * Vectors for use in the iterative solution process
	 */
	private Vector w, u, r;

	/**
	 * Vectors spanning the subspace
	 */
	private Vector[] v;

	/**
	 * Restart vector
	 */
	private DenseVector s;

	/**
	 * Hessenberg matrix
	 */
	private DenseMatrix H;

	/**
	 * Givens rotations for the QR factorization
	 */
	private GivensRotation[] rotation;

	/**
	 * Constructor for GMRES. Uses the given vector as template for creating
	 * scratch vectors. Typically, the solution or the right hand side vector
	 * can be passed, and the template is not modified. The iteration is
	 * restarted every 30 iterations
	 *
	 * @param template
	 *            Vector to use as template for the work vectors needed in the
	 *            solution process
	 */
	public GMRES(Vector template) {
		this(template, 30);
	}

	/**
	 * Constructor for GMRES. Uses the given vector as template for creating
	 * scratch vectors. Typically, the solution or the right hand side vector
	 * can be passed, and the template is not modified
	 *
	 * @param template
	 *            Vector to use as template for the work vectors needed in the
	 *            solution process
	 * @param restart
	 *            GMRES iteration is restarted after this number of iterations
	 */
	public GMRES(Vector template, int restart) {
		w = template.copy();
		u = template.copy();
		r = template.copy();
		setRestart(restart);
	}

	/**
	 * Sets the restart parameter
	 *
	 * @param restart
	 *            GMRES iteration is restarted after this number of iterations
	 */
	public void setRestart(int restart) {
		this.restart = restart;
		if (restart <= 0)
			throw new IllegalArgumentException(
				"restart must be a positive integer");

		s = new DenseVector(restart + 1);
		H = new DenseMatrix(restart + 1, restart);
		rotation = new GivensRotation[restart + 1];

		v = new Vector[restart + 1];
		for (int i = 0; i < v.length; ++i)
			v[i] = r.copy().zero();
	}

	public<T extends VectorMultiplyable> linalg.Vector solve(T preconditioner, VectorMultiplyable A,
	                                                         linalg.Vector b,
	                                                         double tol)
		throws IterativeSolverNotConvergedException
	{
		linalg.Vector x = new linalg.DenseVector(b.getShape().get(0));
		linalg.Vector r_ = preconditioner.mvMul(b.sub(A.mvMul(x)));
		linalg.Vector w;
		double normr = r_.euclidianNorm();
		linalg.Vector v_[] = new linalg.DenseVector[v.length];
		// Outer iteration
		for (iter.setFirst(); normr > tol; iter.next()) {

			v_[0] = r_.mul(1./normr);
			s.zero().set(0, normr);
			int i = 0;
			// Inner iteration
			for (; i < restart && normr>tol; i++, iter
				.next()) {
				w = preconditioner.mvMul(A.mvMul(v_[i]));
				for (int k = 0; k <= i; k++) {
					H.set(k, i, w.inner(v_[k]));
					w = v_[k].mul(-H.get(k, i)).add(w);
				}
				H.set(i + 1, i, w.euclidianNorm());
				v_[i+1] = w.mul(1. / H.get(i + 1, i));
				// QR factorization of H using Givens rotations
				for (int k = 0; k < i; ++k)
					rotation[k].apply(H, i, k, k + 1);

				rotation[i] = new GivensRotation(H.get(i, i), H.get(i + 1, i));
				rotation[i].apply(H, i, i, i + 1);
				rotation[i].apply(s, i, i + 1);
			}

			// Update solution in current subspace
			new UpperTriangDenseMatrix(H, i, false).solve(s, s);
			for (int j = 0; j < i; j++)
				x = x.add(v_[j].mul(s.get(j)));

			r_ = b.sub(A.mvMul(x));
			normr = r_.euclidianNorm();
		}

		return x;
	}

	public linalg.Vector solve(VectorMultiplyable A,
	                                                          linalg.Vector b,
	                                                          double tol)
		throws IterativeSolverNotConvergedException
	{
		linalg.Vector x = new linalg.DenseVector(b.getShape().get(0));
		x.set(Math.random(),0);
		linalg.Vector r_ = b.sub(A.mvMul(x));
		linalg.Vector w;
		double normr = r_.euclidianNorm();
		linalg.Vector v_[] = new linalg.Vector[v.length];
		// Outer iteration
		for (iter.setFirst(); normr > tol; iter.next()) {

			v_[0] = r_.mul(1./normr);
			s.zero().set(0, normr);
			int i = 0;
			// Inner iteration
			for (; i < restart && normr>tol; i++, iter
				.next()) {
				w = A.mvMul(v_[i]);
				for (int k = 0; k <= i; k++) {
					H.set(k, i, w.inner(v_[k]));
					w = v_[k].mul(-H.get(k, i)).add(w);
				}
				H.set(i + 1, i, w.euclidianNorm());
				v_[i+1] = w.mul(1. / H.get(i + 1, i));
				// QR factorization of H using Givens rotations
				for (int k = 0; k < i; ++k)
					rotation[k].apply(H, i, k, k + 1);

				rotation[i] = new GivensRotation(H.get(i, i), H.get(i + 1, i));
				rotation[i].apply(H, i, i, i + 1);
				rotation[i].apply(s, i, i + 1);
			}

			// Update solution in current subspace
			new UpperTriangDenseMatrix(H, i, false).solve(s, s);
			for (int j = 0; j < i; j++)
				x = x.add(v_[j].mul(s.get(j)));

			r_ = b.sub(A.mvMul(x));
			normr = r_.euclidianNorm();
			System.out.println(normr);
		}

		return x;
	}
	@Override
	public Vector solve(Matrix matrix, Vector vector, Vector vector1) throws IterativeSolverNotConvergedException
	{
		return null;
	}
}