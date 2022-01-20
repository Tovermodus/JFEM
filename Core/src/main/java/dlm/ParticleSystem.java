package dlm;

import linalg.Matrix;
import linalg.Vector;

public class ParticleSystem
{
	public final Matrix massMatrix;
	public final Matrix elasticityMatrix;
	public final Matrix semiImplicitMatrix;
	public final Matrix lagrangeMatrix;
	public final Matrix lagrangeBackgroundMatrix;
	public final Vector forceRhs;
	public final Vector accelerationRhs;
	public final Vector lagrangeRhs;
	
	public ParticleSystem(final Matrix massMatrix,
	                      final Matrix elasticityMatrix,
	                      final Matrix semiImplicitMatrix,
	                      final Matrix lagrangeMatrix,
	                      final Matrix lagrangeBackgroundMatrix,
	                      final Vector forceRhs,
	                      final Vector accelerationRhs,
	                      final Vector lagrangeRhs)
	{
		this.massMatrix = massMatrix;
		this.elasticityMatrix = elasticityMatrix;
		this.semiImplicitMatrix = semiImplicitMatrix;
		this.lagrangeMatrix = lagrangeMatrix;
		this.lagrangeBackgroundMatrix = lagrangeBackgroundMatrix;
		this.forceRhs = forceRhs;
		this.accelerationRhs = accelerationRhs;
		this.lagrangeRhs = lagrangeRhs;
	}
}
