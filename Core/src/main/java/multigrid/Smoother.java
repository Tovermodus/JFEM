package multigrid;

import linalg.Vector;
import linalg.VectorMultiplyable;

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
}
