package multigrid;

import linalg.Vector;
import linalg.VectorMultiplyable;

public class MGPReconditioner
	implements VectorMultiplyable
{
	final MGInterface mg;
	
	public MGPReconditioner(final MGInterface mg)
	{
		this.mg = mg;
	}
	
	@Override
	public int getVectorSize()
	{
		return mg.getFinestSystem()
		         .getVectorSize();
	}
	
	@Override
	public int getTVectorSize()
	{
		return mg.getFinestSystem()
		         .getVectorSize();
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
