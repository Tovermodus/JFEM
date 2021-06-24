package mixed;

import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.Tensor;

public class PressureValue extends MixedValue
{
	public PressureValue(double pressure)
	{
		super(3);
		setPressure(pressure);
	}
	
	public PressureValue(double pressure, int domainDimension)
	{
		super(domainDimension);
		setPressure(pressure);
	}
	
	@Override
	public double at(int... coordinates)
	{
		if (coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] != 0)
			throw new IllegalStateException("Not a velocity vector");
		return getPressure();
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		if (coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] != 0)
			throw new IllegalStateException("Not a velocity vector");
		setPressure(value);
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		if (coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] != 0)
			throw new IllegalStateException("Not a velocity vector");
		addPressure(value);
	}
	
	@Override
	public MixedValue add(Tensor other)
	{
		if (other instanceof PressureValue)
			return new PressureValue(getPressure() + ((PressureValue) other).getPressure());
		if (other instanceof VelocityValue)
		{
			MixedValue ret = new MixedValue(((VelocityValue) other).getDomainDimension());
			ret.setPressure(getPressure());
			ret.setVelocity(((VelocityValue) other).getVelocity());
			return ret;
		} else if (other instanceof MixedValue)
		{
			MixedValue ret = new MixedValue(((MixedValue) other).getDomainDimension());
			ret.setPressure(((MixedValue) other).getPressure() + getPressure());
			ret.setVelocity(((MixedValue) other).getVelocity());
			return ret;
		} else
			return super.add(other);
	}
	
	@Override
	public PressureValue mul(double scalar)
	{
		return new PressureValue(scalar);
	}
	
	
	@Override
	public CoordinateVector getVelocity()
	{
		throw new IllegalStateException("Not a Velocity vector");
	}
	
	@Override
	public IntCoordinates getVelocityShape()
	{
		throw new IllegalStateException("Not a Velocity vector");
	}
	
	@Override
	public int getDomainDimension()
	{
		throw new IllegalStateException("Not a Velocity vector");
	}
	
	@Override
	public void setVelocity(CoordinateVector velocity)
	{
		throw new IllegalStateException("Not a Velocity vector");
	}
	
	@Override
	public void addVelocity(CoordinateVector velocity)
	{
		throw new IllegalStateException("Not a Velocity vector");
	}
	
	@Override
	public String toString()
	{
		return "Pressure: " + getPressure();
	}
}