package dlm;

import linalg.Matrix;
import linalg.Vector;

public class FluidSystem
{
	public final Matrix massMatrix;
	public final Matrix flowMatrix;
	public final Matrix semiImplicitMatrix;
	public final Vector forceRhs;
	public final Vector accelerationRhs;
	
	public FluidSystem(final Matrix massMatrix, final Matrix flowMatrix,
	                   final Matrix semiImplicitMatrix,
	                   final Vector forceRhs,
	                   final Vector accelerationRhs)
	{
		this.massMatrix = massMatrix;
		this.flowMatrix = flowMatrix;
		this.semiImplicitMatrix = semiImplicitMatrix;
		this.forceRhs = forceRhs;
		this.accelerationRhs = accelerationRhs;
	}
}
