package dlm;
//GASPAR, NOTAY, OOSTERLEE 2014

import linalg.*;
import multigrid.ForwardBackwardGaussSeidelSmoother;
import multigrid.Smoother;

public class UzawaStokesSmoother3
	implements Smoother
{
	final double omega;
	final int runs;
	final int vSize;
	SparseMatrix A;
	SparseMatrix B;
	SparseMatrix C;
	VectorMultiplyable twoGs = null;
	VectorMultiplyable Sinv;
	
	public UzawaStokesSmoother3(final int runs, final double omega, final int vSize)
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
		
		//return iterate;
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
			B = ((SparseMatrix) Operator).slice(new IntCoordinates(vel_size, 0),
			                                    new IntCoordinates(tot_size, vel_size));
			A = ((SparseMatrix) Operator).slice(new IntCoordinates(0, 0),
			                                    new IntCoordinates(vel_size, vel_size));
			C = ((SparseMatrix) Operator).slice(
				new IntCoordinates(vel_size, vel_size),
				new IntCoordinates(tot_size, tot_size));
			
			final DenseVector AD =
				A.diag();
			SparseMatrix alphaDInv = new SparseMatrix(A.getShape());
			for (int i = 0; i < AD.size(); i++)
			{
				alphaDInv.set(1. / AD.at(i), i, i);
			}
			alphaDInv =
				new SparseMatrix(new BlockDenseMatrix(A,
				                                      A.getCols() / 2).getInvertedDiagonalMatrix());
			final double eigA = VectorMultiplyable.concatenate(alphaDInv, A)
			                                      .powerIterationNonSymmetric();
			System.out.println("AEIG" + eigA);
			final SparseMatrix S = C.add(B.mmMul(alphaDInv)
			                              .mtMul(B)
			                              .mul(eigA));
			final SparseMatrix ADiag = A.mmMul(alphaDInv);
			Sinv = new ForwardBackwardGaussSeidelSmoother(50, S).asPreconditioner(S);
			final VectorMultiplyable ADgs = new ForwardBackwardGaussSeidelSmoother(400,
			                                                                       ADiag).asPreconditioner(
				ADiag);
			final VectorMultiplyable Ags = new ForwardBackwardGaussSeidelSmoother(600,
			                                                                      A).asPreconditioner(
				A);
			final SparseMatrix finalAlphaDInv = alphaDInv;
			twoGs = new VectorMultiplyable()
			{
				@Override
				public int getVectorSize()
				{
					return A.getVectorSize();
				}
				
				@Override
				public int getTVectorSize()
				{
					return A.getTVectorSize();
				}
				
				@Override
				public Vector mvMul(final Vector vector)
				{
					final Vector y = ADgs.mvMul(vector);
					return finalAlphaDInv.mvMul(y);
				}
				
				@Override
				public Vector tvMul(final Vector vector)
				{
					throw new UnsupportedOperationException("not implemented yet");
				}
			};
//				VectorMultiplyable
//				.concatenate(A.getDiagonalMatrix(),
//				             new ForwardBackwardGaussSeidelSmoother(100,
//				                                                    ADiag).asPreconditioner(ADiag));
		}
		final VectorMultiplyable Minv = twoGs;
		for (int i = 0; i < runs; i++)
		{
			final Vector min = Minv.mvMul(f.sub(A.mvMul(u))
			                               .sub(B.tvMul(p)));
			System.out.println(prefix + " MINRES " + A.mvMul(min)
			                                          .sub(f.sub(A.mvMul(u))
			                                                .sub(B.tvMul(p)))
			                                          .euclidianNorm());
			u = u.add(min
				          .mul(1));
			final Vector newp = Sinv.mvMul(B.mvMul(u)
			                                .sub(C.mvMul(p))
			                                .sub(g))
			                        .add(p);
			u = u.sub(Minv.mvMul(B.tvMul(newp.sub(p))));
//			final Vector p = p.add(g.sub(B.mvMul(u))
//			                           .sub(C.mvMul(p))
//			                           .mul(1));
			p = newp;
			if (verbose)
			{
				System.out.println(prefix + "1 " + A.mvMul(u)
				                                    .add(B.tvMul(p))
				                                    .sub(f)
				                                    .euclidianNorm());
				System.out.println(prefix + "2 " + B.mvMul(u)
				                                    .add(C.mvMul(p)
				                                          .sub(g))
				                                    .euclidianNorm());
				System.out.println(prefix + "3 " + Operator.mvMul(DenseVector.concatenate(u, p))
				                                           .sub(rhs)
				                                           .euclidianNorm());
			}
			if (verbose)
				System.out.println(prefix + "iter " + i);
		}
		return DenseVector.concatenate(u, p);
	}
}
