package dlm;
//GASPAR, NOTAY, OOSTERLEE 2014

import linalg.*;
import multigrid.Smoother;

public class UzawaStokesSmoother2
	implements Smoother
{
	final double omega;
	final int runs;
	final int vSize;
	
	VectorMultiplyable twoGs = null;
	
	public UzawaStokesSmoother2(final int runs, final double omega, final int vSize)
	{
		this.vSize = vSize;
		this.omega = omega;
		this.runs = runs;
	}
	
	@Override
	public Vector smooth(final VectorMultiplyable Operator,
	                     final Vector rhs,
	                     final Vector iterate,
	                     final boolean verbose)
	{
		
		//return iterate;
		final int vel_size = vSize;
		final int tot_size = Operator.getVectorSize();
		Vector u = iterate.slice(0, vel_size);
		Vector p = iterate.slice(vel_size, tot_size);
		final Vector f = rhs.slice(0, vel_size);
		final Vector g = rhs.slice(vel_size, tot_size);
		final SparseMatrix B = ((SparseMatrix) Operator).slice(new IntCoordinates(vel_size, 0),
		                                                       new IntCoordinates(tot_size, vel_size));
		final SparseMatrix A = ((SparseMatrix) Operator).slice(new IntCoordinates(0, 0),
		                                                       new IntCoordinates(vel_size, vel_size));
		final SparseMatrix C = ((SparseMatrix) Operator).slice(
			new IntCoordinates(vel_size, vel_size),
			new IntCoordinates(tot_size, tot_size));
		if (verbose)
			System.out.println("    init");
		if (twoGs == null)
			twoGs = new VectorMultiplyable()
			{
				final UpperTriangularSparseMatrix U = new UpperTriangularSparseMatrix(A);
				final LowerTriangularSparseMatrix L = new LowerTriangularSparseMatrix(A);
				final SparseMatrix D = A.getDiagonalMatrix();
				
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
					return L.solve(D.mvMul(U.solve(vector)));
				}
				
				@Override
				public Vector tvMul(final Vector vector)
				{
					throw new UnsupportedOperationException("not implemented yet");
				}
			};
		final VectorMultiplyable Minv = A.inverse();//twoGs;
		for (int i = 0; i < runs; i++)
		{
			u = u.add(Minv.mvMul(f.sub(A.mvMul(u))
			                      .sub(B.tvMul(p))));
			p = B.mvMul(u)
			     .sub(C.mvMul(p))
			     .sub(g)
			     .mul(omega)
			     .add(p);
			if (verbose)
			{
				System.out.println("    1 " + A.mvMul(u)
				                               .add(B.tvMul(p))
				                               .sub(f)
				                               .euclidianNorm());
				System.out.println("    2 " + B.mvMul(u)
				                               .add(C.mvMul(p))
				                               .euclidianNorm());
				System.out.println("    3 " + Operator.mvMul(DenseVector.concatenate(u, p))
				                                      .sub(rhs)
				                                      .euclidianNorm());
			}
			if (verbose)
				System.out.println("    iter " + i);
		}
		return DenseVector.concatenate(u, p);
	}
}
