package mixed;

import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.Tensor;

public class PressureValue
	extends MixedValue
{
	public PressureValue(final double pressure, final int domainDimension)
	{
		super(pressure, new CoordinateVector(domainDimension));
	}
	
	@Override
	public double at(final int... coordinates)
	{
		if (coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] != 0)
			throw new IllegalStateException("Not a velocity vector");
		return getPressure();
	}
	
	@Override
	public void set(final double value, final int... coordinates)
	{
		if (coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] != 0)
			throw new IllegalStateException("Not a velocity vector");
		setPressure(value);
	}
	
	@Override
	public void add(final double value, final int... coordinates)
	{
		if (coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] != 0)
			throw new IllegalStateException("Not a velocity vector");
		addPressure(value);
	}
	
	@Override
	public MixedValue add(final Tensor other)
	{
		if (other instanceof PressureValue)
			return new PressureValue(getPressure() + ((PressureValue) other).getPressure(),
			                         getDomainDimension());
		if (other instanceof VelocityValue)
		{
			final MixedValue ret = new MixedValue(((VelocityValue) other).getDomainDimension());
			ret.setPressure(getPressure());
			ret.setVelocity(((VelocityValue) other).getVelocity());
			return ret;
		} else if (other instanceof MixedValue)
		{
			final MixedValue ret = new MixedValue(((MixedValue) other).getDomainDimension());
			ret.setPressure(((MixedValue) other).getPressure() + getPressure());
			ret.setVelocity(((MixedValue) other).getVelocity());
			return ret;
		} else
			return super.add(other);
	}
	
	@Override
	public PressureValue mul(final double scalar)
	{
		return new PressureValue(getPressure() * scalar,
		                         getDomainDimension());
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
	public void setVelocity(final CoordinateVector velocity)
	{
		throw new IllegalStateException("Not a Velocity vector");
	}
	
	@Override
	public void addVelocity(final CoordinateVector velocity)
	{
		throw new IllegalStateException("Not a Velocity vector");
	}
	
	@Override
	public String toString()
	{
		return "Pressure: " + getPressure();
	}
}
