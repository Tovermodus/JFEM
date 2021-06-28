package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.*;

public class QkQkFunction extends MixedShapeFunction<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction,
	ContinuousTPVectorFunction>
{
	public QkQkFunction(@NotNull ContinuousTPShapeFunction pressureFunction)
	{
		super(pressureFunction);
	}
	
	public QkQkFunction(@NotNull ContinuousTPVectorFunction velocityFunction)
	{
		super(velocityFunction);
	}int globalIndex;
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	public void setGlobalIndex(int index)
	{
		globalIndex = index;
	}
}
