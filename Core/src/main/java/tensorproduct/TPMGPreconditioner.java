package tensorproduct;

import linalg.Vector;
import linalg.VectorMultiplyable;

public class TPMGPreconditioner
	implements VectorMultiplyable
{
	final TPMultiGridSpace<?, ?, ?, ?, ?> mg;
	
	public TPMGPreconditioner(final TPMultiGridSpace<?, ?, ?, ?, ?> mg)
	{
		this.mg = mg;
	}
	
	@Override
	public int getVectorSize()
	{
		return mg.finest_system.getVectorSize();
	}
	
	@Override
	public int getTVectorSize()
	{
		return mg.finest_system.getVectorSize();
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
