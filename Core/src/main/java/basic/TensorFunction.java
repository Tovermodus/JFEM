package basic;

import linalg.DoubleTensor;

public abstract class TensorFunction
{
	private int[] dimensions;
	
	public abstract DoubleTensor value(DoubleTensor pos);

	public abstract DoubleTensor derivative(DoubleTensor pos);
	
	public int[] getDimensions()
	{
		return dimensions;
	}
}
