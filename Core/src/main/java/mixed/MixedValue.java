package mixed;

import linalg.CoordinateVector;
import linalg.Tensor;
import linalg.Vector;

import java.util.List;

public class MixedValue implements Vector
{
	private boolean isVelocity;
	private boolean isPure;
	private double pressure;
	private CoordinateVector velocity;
	private MixedValue()
	{
		isVelocity = false;
		isPure = false;
		pressure = 0;
		velocity = new CoordinateVector(3);
	}
	public static MixedValue pressureValue(double p)
	{
		MixedValue vec = new MixedValue();
		vec.isVelocity = false;
		vec.isPure = true;
		vec.setPressure(p);
		return vec;
	}
	public static MixedValue velocityValue(CoordinateVector c)
	{
		MixedValue vec = new MixedValue();
		vec.isVelocity = true;
		vec.isPure = true;
		vec.setVelocity(new CoordinateVector(c));
		return vec;
	}
	public static MixedValue velocityValue(double ... values)
	{
		MixedValue vec = new MixedValue();
		vec.isVelocity = true;
		vec.isPure = true;
		vec.setVelocity(CoordinateVector.fromValues(values));
		return vec;
	}
	@Override
	public double at(int... coordinates)
	{
		
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if(coordinates[0] == 0)
			return getPressure();
		else
			return getVelocity().at(coordinates[0]-1);
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if(coordinates[0] == 0)
			setPressure(pressure);
		else
			getVelocity().set(value,coordinates[0]-1);
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if(coordinates[0] == 0)
			setPressure(getPressure()+value);
		else
			getVelocity().add(value,coordinates[0]-1);
	}
	
	@Override
	public MixedValue add(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		if(!(other instanceof MixedValue))
			throw new IllegalArgumentException("Cant add MixedValueVector to different Vector");
		MixedValue ret = new MixedValue();
		ret.setPressure(getPressure()+((MixedValue) other).getPressure());
		ret.setVelocity(getVelocity().add(((MixedValue) other).getVelocity()));
		if(isPure() && ((MixedValue) other).isPure())
		{
			ret.isPure = true;
			if(isVelocity() && ((MixedValue) other).isVelocity())
				ret.isVelocity = true;
			if(isPressure() && ((MixedValue) other).isPressure())
				ret.isVelocity = false;
		}
		return ret;
	}
	
	@Override
	public MixedValue mul(double scalar)
	{
		MixedValue ret = new MixedValue();
		ret.isPure = isPure();
		ret.isVelocity = isVelocity();
		ret.setPressure(getPressure() * scalar);
		ret.setVelocity(getVelocity().mul(scalar));
		return ret;
	}
	
	private CoordinateVector getVelocity()
	{
		return velocity;
	}
	
	private double getPressure()
	{
		return pressure;
	}
	
	
	@Override
	public int getSparseEntryCount()
	{
		return getLength();
	}
	
	@Override
	public boolean isSparse()
	{
		return false;
	}
	
	@Override
	public List<Integer> getShape()
	{
		return List.of(1+velocity.getLength());
	}
	public List<Integer> getVelocityShape()
	{
		return getVelocity().getShape();
	}
	public int getVelocityLength()
	{
		return getVelocity().getLength();
	}
	
	public boolean isVelocity()
	{
		return isVelocity;
	}
	public boolean isPressure()
	{
		return !isVelocity&&isPure;
	}
	
	public void setPressure(double pressure)
	{
		this.pressure = pressure;
	}
	
	public void setVelocity(CoordinateVector velocity)
	{
		this.velocity = velocity;
	}
	
	public boolean isPure()
	{
		return isPure;
	}
}
