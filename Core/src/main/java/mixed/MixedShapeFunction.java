package mixed;

import basic.*;
import org.jetbrains.annotations.NotNull;

public abstract class MixedShapeFunction<CT extends Cell<CT,FT>,
	FT extends Face<CT,FT>, ST extends MixedShapeFunction<CT,FT,ST, PF, VF>, PF extends ScalarFunction,
	VF extends VectorFunction> extends MixedFunction<PF,
	VF> implements ShapeFunction<CT,
	FT,ST,MixedValue,MixedGradient,MixedHessian>
{
	protected int globalIndex;
	
	public MixedShapeFunction(@NotNull PF pressureFunction)
	{
		super(pressureFunction);
	}
	
	@Override
	public void setGlobalIndex(int index)
	{
		globalIndex = index;
	}
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
}
