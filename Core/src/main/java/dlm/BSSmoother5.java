package dlm;
//GASPAR, NOTAY, OOSTERLEE 2014

import linalg.*;
import mixed.*;
import multigrid.ForwardBackwardGaussSeidelSmoother;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class BSSmoother5
	implements Smoother
{
	final double alpha;
	final double omega;
	final int runs;
	final int vSize;
	
	VectorMultiplyable Ainv = null;
	SparseMatrix A;
	SparseMatrix B;
	SparseMatrix C;
	SparseMatrix ADiaginv;
	VectorMultiplyable Sinv = null;
	MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian>
		mg;
	
	public BSSmoother5(final int runs, final double alpha, final double omega, final int vSize,
	                   final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian> mg)
	{
		this.vSize = vSize;
		this.alpha = alpha;
		this.omega = omega;
		this.runs = runs;
		this.mg = mg;
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
		if (Sinv == null)
		{
			final SparseMatrix[] blocks = ((SparseMatrix) Operator).partition(new IntCoordinates(vel_size,
			                                                                                     vel_size));
			B = blocks[2];
			A = blocks[0];
			C = blocks[3];
			
			ADiaginv =// new SparseMatrix(A.inverse());
				SparseMatrix.identity(vel_size)
				            .mul(1. / alpha);
			Ainv = new ForwardBackwardGaussSeidelSmoother(5, A).asPreconditioner(A);
			final SparseMatrix S = C.sub(B.mmMul((SparseMatrix) ADiaginv)
			                              .mtMul(B));
			final SparseMatrix Spadded = new SparseMatrix(tot_size, tot_size);
			Spadded.addSmallMatrixInPlaceAt(SparseMatrix.identity(vel_size), 0, 0);
			Spadded.addSmallMatrixInPlaceAt(S, vel_size, vel_size);
			final VectorMultiplyable amg = mg.AMGFromMatrix(Spadded,
			                                                (l, M) -> new ForwardBackwardGaussSeidelSmoother(
				                                                5,
				                                                M), false, 3);
			
			Sinv = new VectorMultiplyable()
			{
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
					final DenseVector v = new DenseVector(tot_size);
					v.addSmallVectorAt(vector, vel_size);
					final Vector sinvV = amg.mvMul(v);
					final Vector redV = sinvV.slice(vel_size, tot_size);
					System.out.println(S.mvMul(redV)
					                    .sub(vector)
					                    .euclidianNorm());
					return redV;
				}
				
				@Override
				public Vector tvMul(final Vector vector)
				{
					throw new UnsupportedOperationException("not implemented yet");
				}
			};
		}
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
			System.out.println();
		}
		for (int i = 0; i < runs; i++)
		{
			final Vector ru = f.sub(A.mvMul(u)
			                         .add(B.tvMul(p)));
			final Vector rp = g.sub(B.mvMul(u)
			                         .add(C.mvMul(p)));
			Vector du = ADiaginv.mvMul(ru);
			final Vector dp = Sinv.mvMul(rp.sub(B.mvMul(du)));
			du = du.sub(ADiaginv.mvMul(B.tvMul(dp)));
			u = u.add(du.mul(omega));
			p = p.add(dp.mul(omega));
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
