package linalg;

import io.vavr.Function2;

public class DirectSchur
	extends ExplicitSchurSolver
{
	Function2<SparseMatrix, Vector, Vector> inverseApplier;
	
	public DirectSchur(final BlockSparseMatrix blockMatrix)
	{
		super(blockMatrix);
	}
	
	@Override
	protected Function2<SparseMatrix, Vector, Vector> solveSchur()
	{
		if (inverseApplier == null)
			inverseApplier = new Function2<SparseMatrix, Vector, Vector>()
			{
				DenseMatrix inverse;
				
				@Override
				public Vector apply(final SparseMatrix v1, final Vector v2)
				{
					if (inverse == null)
						inverse = v1.inverse();
					return inverse.mvMul(v2);
				}
			};
		return inverseApplier;
	}
}
