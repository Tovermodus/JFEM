package dlm;
//GASPAR, NOTAY, OOSTERLEE 2014

import linalg.*;
import multigrid.ForwardBackwardGaussSeidelSmoother;
import multigrid.Smoother;

public class BSSmoother2
	implements Smoother
{
	final double omega;
	final int runs;
	final int vSize;
	
	VectorMultiplyable twoGs = null;
	SparseMatrix A;
	SparseMatrix B;
	SparseMatrix C;
	VectorMultiplyable Sinv = null;
	
	public BSSmoother2(final int runs, final double omega, final int vSize)
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
		if (twoGs == null)
		{
			final SparseMatrix[] blocks = ((SparseMatrix) Operator).partition(new IntCoordinates(vel_size,
			                                                                                     vel_size));
			B = blocks[2];
			A = blocks[0];
			C = blocks[3];
			final SparseMatrix ADiag =
				new SparseMatrix(new BlockDenseMatrix(A, A.getRows() / 100)
					                 .getInvertedDiagonalMatrix());
			
			twoGs = new ForwardBackwardGaussSeidelSmoother(5,
			                                               A).asPreconditioner(
				A,
				v -> v.mul(0));
			
			final double maxEig = VectorMultiplyable.concatenate(ADiag, A)
			                                        .powerIterationNonSymmetric();
			System.out.println("Max Eig estimate " + maxEig + " omega " + maxEig);
			final SparseMatrix S = C.add(B.mmMul(ADiag.mul(maxEig))
			                              .mtMul(B));
			Sinv = new ForwardBackwardGaussSeidelSmoother(5,
			                                              S).asPreconditioner(
				S,
				v -> v.mul(0));
			new VectorMultiplyable()
			{
				final UpperTriangularSparseMatrix U = new UpperTriangularSparseMatrix(S);
				final LowerTriangularSparseMatrix L = new LowerTriangularSparseMatrix(S);
				final SparseMatrix D = S.getDiagonalMatrix();
				
				@Override
				public int getVectorSize()
				{
					return S.getVectorSize();
				}
				
				@Override
				public int getTVectorSize()
				{
					return S.getTVectorSize();
				}
				
				@Override
				public Vector mvMul(final Vector vector)
				{
					return L.solve(D.mvMul(U.solve(vector)));
				}
				
				@Override
				public Vector tvMul(final Vector vector)
				{
					throw new UnsupportedOperationException("not implemented yet");
				}
			};
		}
		final VectorMultiplyable Minv = twoGs;
		for (int i = 0; i < runs; i++)
		{
			u = u.add(Minv.mvMul(f.sub(A.mvMul(u))
			                      .sub(B.tvMul(p))));
			final Vector newp = Sinv.mvMul(B.mvMul(u)
			                                .sub(C.mvMul(p))
			                                .sub(g))
			                        .add(p);
			u = u.sub(Minv.mvMul(B.tvMul(newp.sub(p))));
			p = newp;
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
