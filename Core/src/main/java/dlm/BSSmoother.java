package dlm;
//GASPAR, NOTAY, OOSTERLEE 2014

import linalg.*;
import multigrid.ForwardBackwardGaussSeidelSmoother;
import multigrid.Smoother;

public class BSSmoother
	implements Smoother
{
	final double omega;
	final int runs;
	final int vSize;
	
	//VectorMultiplyable twoGs = null;
	VectorMultiplyable Ainv = null;
	SparseMatrix A;
	SparseMatrix B;
	SparseMatrix C;
	VectorMultiplyable Sinv = null;
	
	public BSSmoother(final int runs, final double omega, final int vSize)
	{
		this.vSize = vSize;
		this.omega = omega;
		this.runs = runs;
	}
	
	@Override
	public Vector smooth(final VectorMultiplyable Operator,
	                     final Vector rhs,
	                     final Vector iterate,
	                     final boolean verbose, final String prefix)
	{
		final int vel_size = vSize;
		final int tot_size = Operator.getVectorSize();
		Vector u = iterate.slice(0, vel_size);
		Vector p = iterate.slice(vel_size, tot_size);
		final Vector f = rhs.slice(0, vel_size);
		final Vector g = rhs.slice(vel_size, tot_size);
		if (verbose)
			System.out.println(prefix + "init");
		if (Ainv == null)
		{
			final SparseMatrix[] blocks = ((SparseMatrix) Operator).partition(new IntCoordinates(vel_size,
			                                                                                     vel_size));
			B = blocks[2];
			A = blocks[0];
			C = blocks[3];
			Ainv = new SparseMatrix(A.getShape());
			final Vector ad = A.diag();
			for (int i = 0; i < ad.getLength(); i++)
				((SparseMatrix) Ainv).set(1. / ad.at(i), i, i);
//			twoGs = new ForwardBackwardGaussSeidelSmoother(5,
//			                                               A).asPreconditioner(
//				A,
//				v -> v.mul(0));

//			double maxEig = VectorMultiplyable.concatenate(Ainv, A)
//			                                  .powerIterationNonSymmetric();
//			System.out.println("Max Eig estimate " + maxEig + " omega " + maxEig);
			Ainv = ((SparseMatrix) Ainv).mul(omega);
//			maxEig = VectorMultiplyable.concatenate(Ainv, A)
//			                           .powerIterationNonSymmetric();
//			System.out.println("Max Eig estimate " + maxEig + " omega " + maxEig);
			final SparseMatrix S = C.sub(B.mmMul(((SparseMatrix) Ainv))
			                              .mtMul(B));
			Sinv =// S.inverseNative();
				new ForwardBackwardGaussSeidelSmoother(20,
				                                       S).asPreconditioner(
					S,
					v -> v.mul(0));
//			new VectorMultiplyable()
//			{
//				final UpperTriangularSparseMatrix U = new UpperTriangularSparseMatrix(S);
//				final LowerTriangularSparseMatrix L = new LowerTriangularSparseMatrix(S);
//				final SparseMatrix D = S.getDiagonalMatrix();
//
//				@Override
//				public int getVectorSize()
//				{
//					return S.getVectorSize();
//				}
//
//				@Override
//				public int getTVectorSize()
//				{
//					return S.getTVectorSize();
//				}
//
//				@Override
//				public Vector mvMul(final Vector vector)
//				{
//					return L.solve(D.mvMul(U.solve(vector)));
//				}
//
//				@Override
//				public Vector tvMul(final Vector vector)
//				{
//					throw new UnsupportedOperationException("not implemented yet");
//				}
//			};
		}
		final VectorMultiplyable Minv = Ainv;
		for (int i = 0; i < runs; i++)
		{
//			u = u.add(Minv.mvMul(f.sub(A.mvMul(u))
//			                      .sub(B.tvMul(p))));
//			final Vector newp = Sinv.mvMul(B.mvMul(u)
//			                                .sub(C.mvMul(p))
//			                                .sub(g))
//			                        .add(p);
//			u = u.sub(Minv.mvMul(B.tvMul(newp.sub(p))));
//			p = newp;
			final Vector ru = f.sub(A.mvMul(u))
			                   .sub(B.tvMul(p));
			final Vector rp = g.sub(B.mvMul(u))
			                   .sub(C.mvMul(p));
			final Vector deltp = Sinv.mvMul(rp.sub(B.mvMul(Ainv.mvMul(ru))));
			final Vector deltu = Ainv.mvMul(ru.sub(B.tvMul(deltp)));
			u = u.add(deltu);
			p = p.add(deltp);
			if (verbose)
			{
				
				System.out.println(prefix + A.mvMul(u)
				                             .add(B.tvMul(p))
				                             .sub(f)
				                             .euclidianNorm());
				System.out.println(prefix + B.mvMul(u)
				                             .add(C.mvMul(p))
				                             .sub(g)
				                             .euclidianNorm());
				System.out.println(prefix + Operator.mvMul(DenseVector.concatenate(u, p))
				                                    .sub(rhs)
				                                    .euclidianNorm());
			}
			if (verbose)
				System.out.println(prefix + "iter " + i);
		}
		return DenseVector.concatenate(u, p);
	}
}
