package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.*;

public class RTMixedFunction extends MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction, RTShapeFunction>
{
	public RTMixedFunction(@NotNull ContinuousTPShapeFunction pressureFunction)
	{
		super(pressureFunction);
	}
	
	public RTMixedFunction(@NotNull RTShapeFunction velocityFunction)
	{
		super(velocityFunction);
	}
}
