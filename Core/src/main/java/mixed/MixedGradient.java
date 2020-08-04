package mixed;

import linalg.*;

import java.util.List;

public class MixedGradient implements Matrix
{
	private boolean isVelocity;
	private boolean isPure;
	private CoordinateVector pressureGradient;
	private CoordinateMatrix velocityGradient;
	private MixedGradient()
	{
		isVelocity = false;
		isPure = false;
		pressureGradient = new CoordinateVector(3);
		velocityGradient = new CoordinateMatrix(3,3);
	}
	@Override
	public double at(int... coordinates)
	{
		if(coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if(coordinates[0] == 0)
			return pressureGradient.at(coordinates[1]);
		else
			return velocityGradient.at(coordinates[0]-1,coordinates[1]);
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		
		if(coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if(coordinates[0] == 0)
			pressureGradient.set(value, coordinates[1]);
		else
			velocityGradient.set(value,coordinates[0]-1,coordinates[1]);
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		if(coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if(coordinates[0] == 0)
			pressureGradient.add(value, coordinates[1]);
		else
			velocityGradient.add(value,coordinates[0]-1,coordinates[1]);
	}
	
	@Override
	public MixedGradient add(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		if(!(other instanceof MixedGradient))
			throw new IllegalArgumentException("Cant add MixedGradient to different Vector");
		MixedGradient ret = new MixedGradient();
		ret.setPressureGradient(getPressureGradient().add(((MixedGradient) other).getPressureGradient()));
		ret.setVelocityGradient(getVelocityGradient().add(((MixedGradient) other).getVelocityGradient()));
		if(isPure() && ((MixedGradient) other).isPure())
		{
			ret.isPure = true;
			if(isVelocity() && ((MixedGradient) other).isVelocity())
				ret.isVelocity = true;
			if(isPressure() && ((MixedGradient) other).isPressure())
				ret.isVelocity = false;
		}
		return ret;
	}
	
	private CoordinateMatrix getVelocityGradient()
	{
		return velocityGradient;
	}
	
	private CoordinateVector getPressureGradient()
	{
		return pressureGradient;
	}
	
	@Override
	public MixedGradient mul(double scalar)
	{
		MixedGradient ret = new MixedGradient();
		ret.isPure = isPure();
		ret.isVelocity = isVelocity();
		ret.setVelocityGradient(getVelocityGradient().mul(scalar));
		ret.setPressureGradient(getPressureGradient().mul(scalar));
		return ret;
	}
	
	@Override
	public List<Vector> unfoldDimension(int dimension)
	{
		return null;
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return getPressureGradient().getLength()+getVelocityGradient().getSparseEntryCount();
	}
	
	@Override
	public boolean isSparse()
	{
		return false;
	}
	
	@Override
	public List<Integer> getShape()
	{
		return List.of(1+getVelocityGradient().getRows(),pressureGradient.getLength());
	}
	
	@Override
	public Matrix transpose()
	{
		throw new UnsupportedOperationException("can't transpose mixed gradient");
	}
	
	@Override
	public Vector mvMul(Vector vector)
	{
		throw new UnsupportedOperationException("can't multiply mixed gradient");
	}
	
	@Override
	public Matrix mmMul(Matrix matrix)
	{
		throw new UnsupportedOperationException("can't multiply mixed gradient");
	}
	
	public boolean isPure()
	{
		return isPure;
	}
	
	public boolean isVelocity()
	{
		return isVelocity;
	}
	public boolean isPressure()
	{
		return !isVelocity()&&isPure();
	}
	
	public void setPressureGradient(CoordinateVector pressureGradient)
	{
		this.pressureGradient = pressureGradient;
	}
	
	public void setVelocityGradient(CoordinateMatrix velocityGradient)
	{
		this.velocityGradient = velocityGradient;
	}
}
