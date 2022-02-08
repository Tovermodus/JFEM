package dlm;
//WABRO DISS

import linalg.*;
import multigrid.ForwardBackwardGaussSeidelSmoother;
import multigrid.Smoother;

public class BSSmoother4
	implements Smoother
{
	final double omega;
	final int runs;
	private final DLMLagrangeAMGSolver.DLMAMG
		amgPreconditionerSpace;
	final int vSize;
	
	VectorMultiplyable Ainv;
	VectorMultiplyable Sinv;
	SparseMatrix C;
	SparseMatrix A;
	SparseMatrix B;
	
	public BSSmoother4(final int runs,
	                   final double omega,
	                   final DLMLagrangeAMGSolver.DLMAMG amgPreconditionerSpace,
	                   final int vSize)
	{
		this.amgPreconditionerSpace = amgPreconditionerSpace;
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
		if (A == null)
		{
			final SparseMatrix[] blocks = ((SparseMatrix) Operator).partition(new IntCoordinates(vel_size,
			                                                                                     vel_size));
			B = blocks[2];
			A = blocks[0];
			C = blocks[3];
			if (verbose)
				System.out.println(prefix + "init");
			final var velocityAMG = amgPreconditionerSpace.getVelocityAMG();
			velocityAMG.verbose = false;
			Ainv = velocityAMG;//Ags.asPreconditioner(A);
			final SparseMatrix Smid =
				new SparseMatrix(new BlockDenseMatrix(A,
				                                      A.getCols()).getInvertedDiagonalMatrix());
			final SparseMatrix S = C.add(B.mmMul(Smid)
			                              .mtMul(B))
			                        .mul(10);
			final VectorMultiplyable Svec = C.addVm(VectorMultiplyable.concatenate(B,
			                                                                       VectorMultiplyable.concatenateTranspose(
				                                                                       Ainv,
//				                                                                       Ags.asPreconditioner(
//					                                                                       A),
				                                                                       B)));
			final Smoother Sgs = new ForwardBackwardGaussSeidelSmoother(5, S);
			Sinv =
				new VectorMultiplyable()
				{
					final IterativeSolver it = new IterativeSolver(true);
					
					@Override
					public int getVectorSize()
					{
						return A.getVectorSize();
					}
					
					@Override
					public int getTVectorSize()
					{
						return A.getVectorSize();
					}
					
					@Override
					public Vector mvMul(final Vector vector)
					{
						
						final Vector sol
							= it
							.solvePGMRES(Svec,
							             Sgs.asPreconditioner(Svec),
							             vector,
							             Math.max(1e-7, 1e-1 * vector.euclidianNorm()));
						System.out.println(prefix + "s iterations" + it.iterations);
//					System.out.println(prefix + "sgs " + vector.sub(Svec.mvMul(vector.add(sol)))
//					                                           .euclidianNorm());
						System.out.println();
						return sol;
//					return Sgs.smooth(S, vector, vector, true, prefix + " sgs ");
					}
					
					@Override
					public Vector tvMul(final Vector vector)
					{
						throw new UnsupportedOperationException("not implemented yet");
					}
				};
		}
		for (int i = 0; i < runs; i++)
		{
			u = u.add(Ainv.mvMul(f.sub(A.mvMul(u))
			                      .sub(B.tvMul(p))));
			final Vector newp = Sinv.mvMul(B.mvMul(u)
			                                .sub(C.mvMul(p))
			                                .sub(g))
			                        .add(p);
			u = u.sub(Ainv.mvMul(B.tvMul(newp.sub(p))));
			
			p = newp;
			if (verbose)
			{
				
				System.out.println(prefix + "first " + A.mvMul(u)
				                                        .add(B.tvMul(p))
				                                        .sub(f)
				                                        .euclidianNorm());
				System.out.println(prefix + "second " + B.mvMul(u)
				                                         .add(C.mvMul(p))
				                                         .sub(g)
				                                         .euclidianNorm());
				System.out.println(prefix + "both " + Operator.mvMul(DenseVector.concatenate(u, p))
				                                              .sub(rhs)
				                                              .euclidianNorm());
			}
			if (verbose)
				System.out.println(prefix + "#####iter " + i);
		}
		System.out.println("BS4RETURN");
		return DenseVector.concatenate(u, p);
	}
}
