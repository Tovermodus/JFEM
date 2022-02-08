package linalg;

import scala.Function1;
import scala.Function2;
import scala.Tuple2;

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
				
				@Override
				public Function1<SparseMatrix, Function1<Vector, Vector>> curried()
				{
					return Function2.super.curried();
				}
				
				@Override
				public Function1<Tuple2<SparseMatrix, Vector>, Vector> tupled()
				{
					return Function2.super.tupled();
				}
			};
		return inverseApplier;
	}
}
