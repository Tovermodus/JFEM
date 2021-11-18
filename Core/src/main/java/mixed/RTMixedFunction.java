package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.RTShapeFunction;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class RTMixedFunction
	extends ComposeMixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction,
	RTShapeFunction>
	implements Comparable<RTMixedFunction>
{
	public RTMixedFunction(@NotNull final ContinuousTPShapeFunction pressureFunction)
	{
		super(pressureFunction);
	}
	
	public RTMixedFunction(@NotNull final RTShapeFunction velocityFunction)
	{
		super(velocityFunction);
	}
	
	int globalIndex;
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	public void setGlobalIndex(final int index)
	{
		globalIndex = index;
	}
	
	@Override
	public int hashCode()
	{
		if (hasPressureFunction())
			return getPressureFunction().hashCode() * 7 + 1;
		if (hasVelocityFunction())
			return getVelocityFunction().hashCode() * 7 + 3;
		return 0;
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (obj instanceof RTMixedFunction)
		{
			if (hasPressureFunction() && ((RTMixedFunction) obj).hasPressureFunction())
				return getPressureFunction().equals(((RTMixedFunction) obj).getPressureFunction());
			if (hasPressureFunction() && ((RTMixedFunction) obj).hasVelocityFunction())
				return false;
			if (hasVelocityFunction() && ((RTMixedFunction) obj).hasPressureFunction())
				return false;
			if (hasVelocityFunction() && ((RTMixedFunction) obj).hasVelocityFunction())
				return getVelocityFunction().equals(((RTMixedFunction) obj).getVelocityFunction());
		}
		return false;
	}
	
	@Override
	public String toString()
	{
		String ret = "|  ";
		if (hasPressureFunction())
			ret += getPressureFunction().toString();
		ret += "  |  ";
		if (hasVelocityFunction())
			ret += getVelocityFunction().toString();
		ret += "  |";
		return ret;
	}
	
	@Override
	public int compareTo(@NotNull final RTMixedFunction o)
	{
		if (hasPressureFunction() && o.hasPressureFunction())
			return getPressureFunction().compareTo(o.getPressureFunction());
		else if (hasPressureFunction() && o.hasVelocityFunction())
			return 1;
		else if (hasVelocityFunction() && o.hasPressureFunction())
			return -1;
		else
			return getVelocityFunction().compareTo(o.getVelocityFunction());
	}
}
