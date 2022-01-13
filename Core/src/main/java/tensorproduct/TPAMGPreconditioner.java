package tensorproduct;

import linalg.Vector;
import linalg.VectorMultiplyable;

public class TPAMGPreconditioner
	implements VectorMultiplyable
{
	final TPAlgebraicMultiGridSpace<?, ?, ?, ?, ?> mg;
	
	public TPAMGPreconditioner(final TPAlgebraicMultiGridSpace<?, ?, ?, ?, ?> mg)
	{
		this.mg = mg;
	}
	
	@Override
	public int getVectorSize()
	{
		return mg.matrix.getVectorSize();
	}
	
	@Override
	public int getTVectorSize()
	{
		return mg.matrix.getVectorSize();
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		return mg.vCycle(vector.mul(0), vector);
	}
	
	@Override
	public Vector tvMul(final Vector vector)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
