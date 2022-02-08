package dlm;
//GASPAR, NOTAY, OOSTERLEE 2014

import linalg.*;
import multigrid.ForwardBackwardGaussSeidelSmoother;
import multigrid.Smoother;

public class BSSmoother3
	implements Smoother
{
	final double alpha;
	final double omega;
	final int runs;
	final int vSize;
	SparseMatrix B;
	SparseMatrix A;
	SparseMatrix BT;
	SparseMatrix C;
	DenseVector AD;
	SparseMatrix alphaDInv;
	VectorMultiplyable Sinv;
	SparseMatrix S;
	SparseMatrix[] blocks;
	
	public BSSmoother3(final int runs, final double alpha, final double omega, final int vSize)
	{
		this.omega = omega;
		this.vSize = vSize;
		this.alpha = alpha;
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
		if (Sinv == null)
		{
			blocks = ((SparseMatrix) Operator).partition(new IntCoordinates(vel_size,
			                                                                vel_size));
			B = blocks[2];
			BT = blocks[1];
			A = blocks[0];
			C = blocks[3];
			AD =
				A.diag();
			alphaDInv = new SparseMatrix(A.getShape());
			for (int i = 0; i < AD.size(); i++)
			{
				alphaDInv.set(1. / AD.at(i), i, i);
			}
			final double eigEst = A.powerIterationNonSymmetric();
			System.out.println("max eig est " + eigEst);
			System.out.println("max eig est " + VectorMultiplyable.concatenate(alphaDInv, A)
			                                                      .powerIterationNonSymmetric());
			S = C.sub(B.mmMul(alphaDInv)
			           .mmMul(BT));
			Sinv = new VectorMultiplyable()
			{
				final Smoother smooth = new ForwardBackwardGaussSeidelSmoother(5, S);
				
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
					return smooth.smooth(S, vector, vector, true, prefix + " gs: ");
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
			
			System.out.println(prefix + "init");
			System.out.println(prefix + String.format("%6.3e", A.mvMul(u)
			                                                    .add(BT.mvMul(p))
			                                                    .sub(f)
			                                                    .euclidianNorm()));
			System.out.println(prefix + String.format("%6.3e", B.mvMul(u)
			                                                    .add(C.mvMul(p))
			                                                    .sub(g)
			                                                    .euclidianNorm()));
			System.out.println(prefix + String.format("%6.3e", Operator.mvMul(DenseVector.concatenate(u, p))
			                                                           .sub(rhs)
			                                                           .euclidianNorm()));
		}
		
		for (int i = 0; i < runs; i++)
		{
			final Vector r = A.mvMul(u)
			                  .add(B.tvMul(p))
			                  .sub(f);
			//System.out.println("r" + r.absMaxElement());
			final Vector s = B.mvMul(u)
			                  .add(C.mvMul(p))
			                  .sub(g);
			//System.out.println("s" + s.absMaxElement());
			final Vector phatRhs = s.sub(B.mvMul(alphaDInv.mvMul(r)));
			//System.out.println("pr" + phatRhs.absMaxElement());
			final Vector pdelta = Sinv.mvMul(phatRhs);
			//System.out.println("pd" + pdelta.absMaxElement());
			final Vector udelta = alphaDInv.mvMul(r.sub(B.tvMul(pdelta)));
			//System.out.println("ud" + udelta.absMaxElement());
			p = p.sub(pdelta.mul(omega));
			u = u.sub(udelta.mul(omega));
			if (verbose)
			{
				
				System.out.println(prefix + "iter " + i);
				System.out.println(prefix + String.format("%6.3e", A.mvMul(u)
				                                                    .add(BT.mvMul(p))
				                                                    .sub(f)
				                                                    .euclidianNorm()));
				System.out.println(prefix + String.format("%6.3e", B.mvMul(u)
				                                                    .add(C.mvMul(p))
				                                                    .sub(g)
				                                                    .euclidianNorm()));
				System.out.println(prefix + String.format("%6.3e",
				                                          Operator.mvMul(DenseVector.concatenate(u, p))
				                                                  .sub(rhs)
				                                                  .euclidianNorm()));
			}
		}
		return DenseVector.concatenate(u, p);
	}
}
