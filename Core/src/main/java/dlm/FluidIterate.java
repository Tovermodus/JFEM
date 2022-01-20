package dlm;

import linalg.DenseVector;
import linalg.Vector;

public class FluidIterate
{
	public final Vector current;
	
	public FluidIterate(final Vector currentVelocity,
	                    final Vector currentPressure)
	{
		current = DenseVector.concatenate(currentVelocity, currentPressure);
	}
	
	public FluidIterate(final Vector current)
	{
		this.current = current;
	}
}
