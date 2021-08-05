package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class RTMixedFunction extends MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction,
	RTShapeFunction>
{
	public RTMixedFunction(@NotNull ContinuousTPShapeFunction pressureFunction)
	{
		super(pressureFunction);
	}
	
	public RTMixedFunction(@NotNull RTShapeFunction velocityFunction)
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
}
