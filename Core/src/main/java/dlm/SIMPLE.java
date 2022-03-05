package dlm;
//GASPAR, NOTAY, OOSTERLEE 2014

import linalg.*;
import multigrid.Smoother;

public class SIMPLE
	implements Smoother
{
	final double omega;
	final int runs;
	final int vSize;
	
	VectorMultiplyable Ainv = null;
	VectorMultiplyable ADiaginv = null;
	SparseMatrix A;
	SparseMatrix B;
	SparseMatrix C;
	VectorMultiplyable Sinv = null;
	
	public SIMPLE(final int runs, final double omega, final int vSize)
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
			Ainv = A.inverse();
			ADiaginv = A.getDiagonalMatrix()
			            .inverse();
			final SparseMatrix S = new SparseMatrix(C.sub(B.mmMul(((Matrix) ADiaginv).mtMul(B))));
			Sinv = S.inverse();
		}
		for (int i = 0; i < runs; i++)
		{
			final Vector ru = f.sub(A.mvMul(u)
			                         .add(B.tvMul(p)));
			final Vector rp = g.sub(B.mvMul(u)
			                         .add(C.mvMul(p)));
			Vector du = Ainv.mvMul(ru);
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
