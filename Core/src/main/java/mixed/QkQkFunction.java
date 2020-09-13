package mixed;

import org.jetbrains.annotations.NotNull;
import tensorproduct.*;

public class QkQkFunction extends MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction, ContinuousTPVectorFunction>
{
	public QkQkFunction(@NotNull ContinuousTPShapeFunction pressureFunction)
	{
		super(pressureFunction);
	}
	
	public QkQkFunction(@NotNull ContinuousTPVectorFunction velocityFunction)
	{
		super(velocityFunction);
	}
}
