package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class QkQkFunction
	extends ComposeMixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction,
	ContinuousTPVectorFunction>
	implements Comparable<QkQkFunction>
{
	public QkQkFunction(@NotNull final ContinuousTPShapeFunction pressureFunction)
	{
		super(pressureFunction);
	}
	
	public QkQkFunction(@NotNull final ContinuousTPVectorFunction velocityFunction)
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
		
		if (hasPressureFunction())
			getPressureFunction().setGlobalIndex(index);
		if (hasVelocityFunction())
			getVelocityFunction().setGlobalIndex(index);
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
		if (obj instanceof QkQkFunction)
		{
			if (hasPressureFunction() && ((QkQkFunction) obj).hasPressureFunction())
				return getPressureFunction().equals(((QkQkFunction) obj).getPressureFunction());
			if (hasPressureFunction() && ((QkQkFunction) obj).hasVelocityFunction())
				return false;
			if (hasVelocityFunction() && ((QkQkFunction) obj).hasPressureFunction())
				return false;
			if (hasVelocityFunction() && ((QkQkFunction) obj).hasVelocityFunction())
				return getVelocityFunction().equals(((QkQkFunction) obj).getVelocityFunction());
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
	public int compareTo(@NotNull final QkQkFunction o)
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
