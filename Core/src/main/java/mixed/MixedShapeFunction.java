package mixed;

import basic.Face;
import basic.ShapeFunction;
import basic.Cell;

public abstract class MixedShapeFunction<CT extends Cell<CT,FT,ST>,
	FT extends Face<CT,FT,ST>, ST extends MixedShapeFunction<CT,FT,ST>> extends MixedFunction implements ShapeFunction<CT,
	FT,ST,MixedValue,MixedGradient,MixedHessian>
{
	protected int globalIndex;
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
