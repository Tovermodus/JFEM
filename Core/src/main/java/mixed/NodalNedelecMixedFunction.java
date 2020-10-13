package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.*;

public class NodalNedelecMixedFunction extends MixedShapeFunction<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction,
	NodalNedelecShapeFunction>
{
	public NodalNedelecMixedFunction(@NotNull ContinuousTPShapeFunction pressureFunction)
	{
		super(pressureFunction);
	}
	
	public NodalNedelecMixedFunction(@NotNull NodalNedelecShapeFunction velocityFunction)
	{
		super(velocityFunction);
	}
}
