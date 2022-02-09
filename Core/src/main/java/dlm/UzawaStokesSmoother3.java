package dlm;
//GASPAR, NOTAY, OOSTERLEE 2014

import linalg.*;
import mixed.*;
import multigrid.ForwardBackwardGaussSeidelSmoother;
import multigrid.MGPreconditionerSpace;
import multigrid.RichardsonSmoother;
import multigrid.Smoother;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class UzawaStokesSmoother3
	implements Smoother
{
	final double omega;
	final int runs;
	private final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian>
		mg;
	final int vSize;
	SparseMatrix A;
	SparseMatrix B;
	SparseMatrix C;
	VectorMultiplyable S;
	VectorMultiplyable Ainv = null;
	VectorMultiplyable Sinner;
	VectorMultiplyable Sinv;
	SparseMatrix matrixS;
	RichardsonSmoother Ssmooth;
	
	public UzawaStokesSmoother3(final int runs,
	                            final double omega,
	                            final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian> mg,
	                            final int vSize)
	{
		this.mg = mg;
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
		if (Ainv == null)
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
				                                      A.getCols() / 4).getInvertedDiagonalMatrix());
			final double eigA = VectorMultiplyable.concatenate(alphaDInv, A)
			                                      .powerIterationNonSymmetric();
			System.out.println("AEIG" + eigA);
			
			final SparseMatrix Apadded = new SparseMatrix(((SparseMatrix) Operator).getShape());
			Apadded.addSmallMatrixInPlaceAt(A, 0, 0);
			Apadded.addSmallMatrixInPlaceAt(SparseMatrix.identity(C.getCols()), vel_size, vel_size);
			final var Aamg
				= mg.AMGFromMatrix(Apadded,
				                   (l, m) ->
					                   new ForwardBackwardGaussSeidelSmoother(50, m));
			Ainv = new VectorMultiplyable()
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
					final DenseVector padded = new DenseVector(Operator.getVectorSize());
					padded.addSmallVectorAt(vector, 0);
					Vector defect = padded;
					Vector sol = Aamg.mvMul(defect);
					defect = padded.sub(Apadded.mvMul(sol));
					//System.out.println(prefix + " SGSAAMG  " + defect.euclidianNorm());
					sol = sol.add(Aamg.mvMul(defect));
					sol = sol.slice(0, vel_size);
//					System.out.println(prefix + " SGSAAMG--  " + A.mvMul(sol)
//					                                              .sub(vector)
//					                                              .euclidianNorm());
					return sol;
				}
				
				@Override
				public Vector tvMul(final Vector vector)
				{
					throw new UnsupportedOperationException("not implemented yet");
				}
			};
			Sinner = C.addVm(VectorMultiplyable.concatenate(B,
			                                                VectorMultiplyable.concatenateTranspose(
				                                                Ainv, B)));
			final SparseMatrix padding = new SparseMatrix(C.getCols(), Operator.getVectorSize());
			padding.addSmallMatrixInPlaceAt(SparseMatrix.identity(C.getCols()), 0,
			                                Operator.getVectorSize() - C.getCols());
			final SparseMatrix Scomplement = new SparseMatrix(((SparseMatrix) Operator).getShape());
			Scomplement.addSmallMatrixInPlaceAt(SparseMatrix.identity(A.getCols()), 0, 0);
			S = VectorMultiplyable.transposeConcatenate(padding,
			                                            VectorMultiplyable.concatenate(Sinner, padding))
			                      .addVm(Scomplement);
			
			matrixS = new SparseMatrix(
				//	padding.tmMul(
				C.add(B.mmMul(alphaDInv.mul(eigA))
				       .mtMul(B))
				//                                        .mmMul(padding)
			);
			//	                                  .add(Scomplement));
			final Smoother Sgs =
				new ForwardBackwardGaussSeidelSmoother(2, matrixS);
			final VectorMultiplyable Samg = mg.AMGFromVectorMultipliable(S,
			                                                             SparseMatrix.identity(mg.systems.get(
				                                                                                     0)
			                                                                                             .getVectorSize()),
			                                                             (l, m) ->
			                                                             {
				                                                             if (l == mg.maxLevel())
					                                                             return Sgs;
				                                                             return new RichardsonSmoother(
					                                                             1,
					                                                             0);
			                                                             },
			                                                             true);
			
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
					
					final DenseVector padded = new DenseVector(Operator.getVectorSize());
					padded.addSmallVectorAt(vector, vel_size);
					return Samg.mvMul(padded)
					           .slice(vel_size, tot_size);
				}
				
				@Override
				public Vector tvMul(final Vector vector)
				{
					throw new UnsupportedOperationException("not implemented yet");
				}
			};
			Ssmooth = new RichardsonSmoother(1, 3,
			                                 new ForwardBackwardGaussSeidelSmoother(30,
			                                                                        matrixS).asPreconditioner(
				                                 matrixS));
		}
		for (int i = 0; i < runs; i++)
		{
			final Vector min = Ainv.mvMul(f.sub(A.mvMul(u))
			                               .sub(B.tvMul(p)));
			System.out.println(prefix + " AMINRES " + A.mvMul(min)
			                                           .sub(f.sub(A.mvMul(u))
			                                                 .sub(B.tvMul(p)))
			                                           .euclidianNorm());
			u = u.add(min.mul(1));
			final Vector newp = p.add(Ssmooth.smooth(Sinner, B.mvMul(u)
			                                                  .sub(C.mvMul(p))
			                                                  .sub(g), p.mul(0), true, "RICHH "));
//				Sinv.mvMul(B.mvMul(u)
//			                                .sub(C.mvMul(p))
//			                                .sub(g));
			System.out.println(prefix + " SMINRES " + Sinner.mvMul(newp.sub(p))
			                                                .sub(B.mvMul(u)
			                                                      .sub(C.mvMul(p))
			                                                      .sub(g))
			                                                .euclidianNorm());
			u = u.sub(Ainv.mvMul(B.tvMul(newp.sub(p))));
			//p = p.add(newp);
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
