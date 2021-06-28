package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.*;

public class MixedNedelecFunction extends MixedShapeFunction<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction,
	NedelecShapeFunction>
{
	public MixedNedelecFunction(@NotNull ContinuousTPShapeFunction pressureFunction)
	{
		super(pressureFunction);
	}
	
	public MixedNedelecFunction(@NotNull NedelecShapeFunction velocityFunction)
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
