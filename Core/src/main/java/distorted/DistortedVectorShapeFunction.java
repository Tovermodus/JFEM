package distorted;

import basic.SingleComponentVectorShapeFunction;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import org.jetbrains.annotations.NotNull;

public class DistortedVectorShapeFunction extends SingleComponentVectorShapeFunction<DistortedCell, DistortedFace, DistortedShapeFunction>
	implements Comparable<DistortedVectorShapeFunction>

{
	public DistortedVectorShapeFunction(final DistortedCell supportCell, final int polynomialDegree, final int localIndex)
	{
		super(supportCell, polynomialDegree, localIndex, DistortedShapeFunction.class);
	}
	
	public DistortedVectorShapeFunction(final DistortedShapeFunction componentFunction, final int component)
	{
		super(componentFunction, component);
	}
	
	@Override
	public int compareTo(@NotNull final DistortedVectorShapeFunction other)
	{
		if (this.getComponent() != other.getComponent())
			return Integer.compare(this.getComponent(), other.getComponent());
		return this.getComponentFunction().compareTo(other.getComponentFunction());
	}
}
