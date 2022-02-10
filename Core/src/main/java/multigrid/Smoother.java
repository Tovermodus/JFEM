package multigrid;

import linalg.Vector;
import linalg.VectorMultiplyable;

import java.util.function.Function;

public interface Smoother
{
	default Vector smooth(final VectorMultiplyable Operator, final Vector rhs, final Vector iterate)
	{
		return smooth(Operator, rhs, iterate, false, "");
	}
	
	Vector smooth(VectorMultiplyable Operator, Vector rhs, Vector iterate, boolean verbose, String prefix);
	
	default VectorMultiplyable asPreconditioner(final VectorMultiplyable op)
	{
		return new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return op.getVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return op.getTVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				final Vector defect = vector.sub(op.mvMul(vector));
				final Vector sol = smooth(op, defect, vector.mul(0));
				return vector.add(sol);
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				throw new UnsupportedOperationException("not implemented yet");
			}
		};
	}
	
	default VectorMultiplyable asPreconditioner(final VectorMultiplyable op,
	                                            final Function<Vector, Vector> rhsToInitial)
	{
		return new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return op.getVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return op.getTVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				final Vector initial = rhsToInitial.apply(vector);
				final Vector defect = vector.sub(op.mvMul(initial));
				final Vector sol = smooth(op, defect, vector.mul(0));
				System.out.println("r" + op.mvMul(initial.add(sol))
				                           .sub(vector)
				                           .euclidianNorm());
				return initial.add(sol);
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				throw new UnsupportedOperationException("not implemented yet");
			}
		};
	}
}
