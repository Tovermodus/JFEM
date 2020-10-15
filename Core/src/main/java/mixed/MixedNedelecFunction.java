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
	}
}
