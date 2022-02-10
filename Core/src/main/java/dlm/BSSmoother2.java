package dlm;
//GASPAR, NOTAY, OOSTERLEE 2014

import linalg.*;
import multigrid.Smoother;

public class BSSmoother2
	implements Smoother
{
	final double omega;
	final int runs;
	final int vSize;
	
	VectorMultiplyable twoGs = null;
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
		//return iterate;
		final int vel_size = vSize;
		final int tot_size = Operator.getVectorSize();
		Vector u = iterate.slice(0, vel_size);
		Vector p = iterate.slice(vel_size, tot_size);
		final Vector f = rhs.slice(0, vel_size);
		final Vector g = rhs.slice(vel_size, tot_size);
		final SparseMatrix O = (SparseMatrix) Operator;
		System.out.println("SYSMMMM" + O.sub(O.transpose())
		                                .absMaxElement());
		final SparseMatrix[] blocks = ((SparseMatrix) Operator).partition(new IntCoordinates(vel_size,
		                                                                                     vel_size));
		final SparseMatrix B = blocks[2];
		final SparseMatrix BT = blocks[1];
		final SparseMatrix A = blocks[0];
		final SparseMatrix C = blocks[3];
//		final SparseMatrix B = ((SparseMatrix) Operator).slice(new IntCoordinates(vel_size, 0),
//		                                                       new IntCoordinates(tot_size, vel_size));
//		final SparseMatrix BT = ((SparseMatrix) Operator).slice(new IntCoordinates(0, vel_size),
//		                                                        new IntCoordinates(vel_size, tot_size));
//		final SparseMatrix A = ((SparseMatrix) Operator).slice(new IntCoordinates(0, 0),
//		                                                       new IntCoordinates(vel_size, vel_size));
//		final SparseMatrix C = ((SparseMatrix) Operator).slice(
//			new IntCoordinates(vel_size, vel_size),
//			new IntCoordinates(tot_size, tot_size));
		if (verbose)
			System.out.println(prefix + "init");
		if (twoGs == null)
			
			if (Sinv == null)
			{
				final SparseMatrix ADiag =
					new SparseMatrix(new BlockDenseMatrix(A, A.getRows() / 100)
						                 .getInvertedDiagonalMatrix());
				
				twoGs =
					//new BlockDenseMatrix(A, A.getRows() / 100).getInvertedDiagonalMatrix();//A.inverse();
					new VectorMultiplyable()
					{
						final SparseMatrix AdA = A.add(A.transpose())
						                          .mul(1. / 2)
						                          .mmMul(ADiag);
						final UpperTriangularSparseMatrix U = new UpperTriangularSparseMatrix(
							AdA);
						final LowerTriangularSparseMatrix L = new LowerTriangularSparseMatrix(
							AdA);
						final SparseMatrix D = AdA.getDiagonalMatrix();
						
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
							final Vector y = L.solve(D.mvMul(U.solve(vector)));
							return ADiag.mvMul(y);
						}
						
						@Override
						public Vector tvMul(final Vector vector)
						{
							throw new UnsupportedOperationException("not implemented yet");
						}
					};
				final double maxEig = VectorMultiplyable.concatenate(ADiag, A)
				                                        .powerIterationNonSymmetric();
				final double omega = -Math.max(1, maxEig * 1.5);
				System.out.println(A.sub(A.transpose())
				                    .absMaxElement());
				System.out.println("Max Eig estimate " + maxEig + " omega " + omega);
				ADiag.mul(1. / omega);
				Sinv = new VectorMultiplyable()
				{
					final SparseMatrix S = C.add(B.mmMul(ADiag)//.mmMul(A.inverse())
					                              .mtMul(B));
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
//			u = u.add(Minv.mvMul(f.sub(A.mvMul(u))
//			                      .sub(B.tvMul(p))));
			final Vector newp = Sinv.mvMul(B.mvMul(u)
			                                .sub(C.mvMul(p))
			                                .sub(g))
			                        .add(p);
			u = u.sub(Minv.mvMul(B.tvMul(newp.sub(p))));
			p = newp;
			if (verbose)
			{
				
				System.out.println(prefix + A.mvMul(u)
				                             .add(BT.mvMul(p))
				                             .sub(f)
				                             .euclidianNorm());
				System.out.println(prefix + B.mvMul(u)
				                             .add(C.mvMul(p))
				                             .euclidianNorm());
				System.out.println(prefix + g.euclidianNorm());
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
