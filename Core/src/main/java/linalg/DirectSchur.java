package linalg;

import scala.Function1;
import scala.Function2;
import scala.Tuple2;

public class DirectSchur
	extends SchurSolver
{
	Function2<DenseMatrix, Vector, Vector> inverseApplier;
	
	public DirectSchur(final BlockSparseMatrix blockMatrix)
	{
		super(blockMatrix);
	}
	
	@Override
	Function2<DenseMatrix, Vector, Vector> solveSchur()
	{
		if (inverseApplier == null)
			inverseApplier = new Function2<DenseMatrix, Vector, Vector>()
			{
				DenseMatrix inverse;
				
				@Override
				public Vector apply(final DenseMatrix v1, final Vector v2)
				{
					if (inverse == null)
						inverse = v1.inverse();
					return inverse.mvMul(v2);
				}
				
				@Override
				public Function1<DenseMatrix, Function1<Vector, Vector>> curried()
				{
					return Function2.super.curried();
				}
				
				@Override
				public Function1<Tuple2<DenseMatrix, Vector>, Vector> tupled()
				{
					return Function2.super.tupled();
				}
			};
		return inverseApplier;
	}
}
