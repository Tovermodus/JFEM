package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.*;

public class NedelecMixedFunction extends MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction, NedelecShapeFunction>
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
