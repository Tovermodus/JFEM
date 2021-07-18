package systems;

import linalg.IntCoordinates;
import linalg.Tensor;

import java.util.List;

public class SystemHessian implements Tensor
{
	@Override
	public double at(int... coordinates)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Tensor add(Tensor other)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Tensor mul(double scalar)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public List<? extends Tensor> unfoldDimension(int dimension)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Tensor slice(IntCoordinates start, IntCoordinates end)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int getSparseEntryCount()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public boolean isSparse()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public IntCoordinates getShape()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public String printFormatted(double... tol)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
