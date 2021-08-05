package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class QkQkFunction extends MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction,
	ContinuousTPVectorFunction> implements Comparable<QkQkFunction>
{
	public QkQkFunction(@NotNull ContinuousTPShapeFunction pressureFunction)
	{
		super(pressureFunction);
	}
	
	public QkQkFunction(@NotNull ContinuousTPVectorFunction velocityFunction)
	{
		super(velocityFunction);
	}
	int globalIndex;
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	public void setGlobalIndex(int index)
	{
		globalIndex = index;
	}
	
	@Override
	public int hashCode()
	{
		if(hasPressureFunction())
			return getPressureShapeFunction().hashCode()*7+1;
		if(hasVelocityFunction())
			return getVelocityShapeFunction().hashCode()*7+3;
		return 0;
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(obj instanceof QkQkFunction)
		{
			if(hasPressureFunction() && ((QkQkFunction) obj).hasPressureFunction())
				return getPressureShapeFunction().equals(((QkQkFunction) obj).getPressureShapeFunction());
			if(hasPressureFunction() && ((QkQkFunction) obj).hasVelocityFunction())
				return false;
			if(hasVelocityFunction() && ((QkQkFunction) obj).hasPressureFunction())
				return false;
			if(hasVelocityFunction() && ((QkQkFunction) obj).hasVelocityFunction())
				return getVelocityShapeFunction().equals(((QkQkFunction) obj).getVelocityShapeFunction());
		}
		return false;
	}
	
	@Override
	public String toString()
	{
		String ret = "|  ";
		if(hasPressureFunction())
			ret += getPressureFunction().toString();
		ret += "  |  ";
		if(hasVelocityFunction())
			ret += getVelocityFunction().toString();
		ret += "  |";
		return ret;
	}
	
	@Override
	public int compareTo(@NotNull QkQkFunction o)
	{
		if(hasPressureFunction() && o.hasPressureFunction())
			return getPressureShapeFunction().compareTo(o.getPressureShapeFunction());
		else if(hasPressureFunction() && o.hasVelocityFunction())
			return 1;
		else if(hasVelocityFunction() && o.hasPressureFunction())
			return -1;
		else
			return getVelocityShapeFunction().compareTo(o.getVelocityShapeFunction());
	}
}
