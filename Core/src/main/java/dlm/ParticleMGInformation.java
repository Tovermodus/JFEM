package dlm;

import linalg.DenseMatrix;

public class ParticleMGInformation
{
	Particle p;
	ParticleIterate iterate;
	DenseMatrix systemInverse;
	
	public ParticleMGInformation(final Particle p)
	{
		this.p = p;
	}
}
