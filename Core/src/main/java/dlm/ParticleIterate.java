package dlm;

import linalg.DenseVector;
import linalg.Vector;

public class ParticleIterate
{
	public final Vector current;
	public final Vector currentLagrange;
	public final Vector last;
	public final Vector lastLagrange;
	
	public ParticleIterate(final Vector currentIterate,
	                       final Vector currentLagrangeIterate,
	                       final Vector lastIterate, final Vector lastLagrangeIterate)
	{
		this.current = currentIterate;
		this.currentLagrange = currentLagrangeIterate;
		this.last = lastIterate;
		this.lastLagrange = lastLagrangeIterate;
	}
	
	public ParticleIterate(final ParticleIterate iterate, final DenseVector newIterate, final DenseVector newLagrangeIterate)
	{
		this.current = newIterate;
		this.currentLagrange = newLagrangeIterate;
		this.last = iterate.current;
		this.lastLagrange = iterate.currentLagrange;
	}
}
