package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.*;

public class NedelecMixedFunction extends MixedShapeFunction<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction,
	NedelecShapeFunction>
{
	public NedelecMixedFunction(@NotNull ContinuousTPShapeFunction pressureFunction)
	{
		super(pressureFunction);
	}
	
	public NedelecMixedFunction(@NotNull NedelecShapeFunction velocityFunction)
	{
		super(velocityFunction);
	}
}
